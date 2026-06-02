module precision 
    ! Double precision kind parameter
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 300)
end module precision

module sim_params
    ! Simulation parameters
    use precision
    implicit none

    type :: params_t
        ! Geometry
        real(kind=dp) :: L_outer, L_inner
        real(kind=dp) :: R_outer, R_inner

        ! Boundary conditions
        real(kind=dp) :: phi_inner  ! Potential on inner conductor

        ! Solver settings
        real(kind=dp) :: h          ! Grid spacing
        real(kind=dp) :: omega      ! SOR relaxation factor
        real(kind=dp) :: tol        ! Convergence tolerance
        integer       :: max_iter   ! Maximum iterations
    end type params_t

    contains 

        function default_params() result(p)
            ! Returns default simulation parameters
            type(params_t) :: p
            p%L_outer = 20.0_dp
            p%L_inner = 5.0_dp
            p%R_outer = 10.0_dp
            p%R_inner = 1.0_dp
            p%phi_inner = 1000.0_dp
        
            p%h        = 0.02_dp
            p%omega    = 1.9_dp
            p%tol      = 1.0e-4_dp
            p%max_iter = 200000
        end function default_params

end module sim_params

module potentials
    use precision 
    use sim_params
    implicit none

    contains 

        ! Analytical solution for coaxial capacitor at z=0
        pure function coax_analytic(r, p) result(val)
        ! Assumes r > R inner (Caller must ensure this)
            real(kind=dp), intent(in) :: r
            type(params_t), intent(in) :: p
            real(kind=dp) :: val 

            ! Calculate potential
            val = (log(p%R_outer) - log(r)) / & 
            (log(p%R_outer) - log(p%R_inner)) * p%phi_inner
        
        end function coax_analytic

        ! SOR update at r=0 (axis)
        pure function update_axis(phi_old, phi_right, phi_up, phi_down, omega) result(phi_new)
            ! Update potential at r=0 using SOR
            real(kind=dp), intent(in) :: phi_old    ! phi(i,0)
            real(kind=dp), intent(in) :: phi_right  ! phi(i,1)
            real(kind=dp), intent(in) :: phi_up     ! phi(i+1,0)
            real(kind=dp), intent(in) :: phi_down   ! phi(i-1,0)

            real(kind=dp), intent(in) :: omega
            real(kind=dp) :: phi_new, U

            ! Compute the update value U
            U = (2.0_dp/3.0_dp) * phi_right + &
            (1.0_dp/6.0_dp) * (phi_up + phi_down)
            
            ! SOR formula
            phi_new = phi_old + omega * (U - phi_old)

        end function update_axis

        ! SOR update at r>0
        pure function update_bulk(phi_old, p_left, p_right, p_up, p_down, inv_r_factor, omega) result(phi_new)

            real(kind=dp), intent(in) :: phi_old        ! phi(i,0)
            real(kind=dp), intent(in) :: p_left         ! phi(i,-1)
            real(kind=dp), intent(in) :: p_right        ! phi(i,1)
            real(kind=dp), intent(in) :: p_up           ! phi(i+1,0)
            real(kind=dp), intent(in) :: p_down         ! phi(i-1,0)

            real(kind=dp), intent(in) :: inv_r_factor   ! Avoid slow division (1/8r)

            real(kind=dp), intent(in) :: omega          
            real(kind=dp) :: phi_new, U                 

            ! Compute the update value U
            U = 0.25_dp * (p_left + p_right + p_up + p_down) + &
            inv_r_factor * (p_right - p_left)
            
            ! SOR formula
            phi_new = phi_old + omega * (U - phi_old)

        end function update_bulk

end module potentials

module grid_data 
    use precision 
    use sim_params
    use potentials
    use mpi 
    implicit none 

    ! Main potential array (allocates with halo cells)
    real(kind=dp), allocatable :: phi(:,:)

    ! Grid dimensions (global)
    integer :: nz, nr                   ! Number of z and r points
    integer :: nz_in, nr_in             ! Number of z and r points for inner conductor

    ! MPI dimensions
    integer :: nr_local                 ! Number of columns (r-indices) for this rank
    integer :: nr_start                 ! Starting global r-index for this rank
    integer :: rank, nprocs             ! MPI rank and number of processes
    integer :: comm, ierr               ! MPI communicator and error code
    integer :: left_rank, right_rank    ! Neighbour ranks

    contains 

        subroutine init_grid(p, mpi_comm)
            ! Initialise grid dimensions and MPI params
            type(params_t), intent(in) :: p
            integer, intent(in) :: mpi_comm
            integer :: i, k_local, k_global, remainder
            real(kind=dp) :: r_val

            comm = mpi_comm
            call MPI_Comm_rank(comm, rank, ierr)
            call MPI_Comm_size(comm, nprocs, ierr)

            ! Global grid size
            nz = nint(p%L_outer / p%h)
            nr = nint(p%R_outer / p%h)
            nz_in = nint(p%L_inner / p%h)
            nr_in = nint(p%R_inner / p%h)

            ! Domain decomposition (split R)
            nr_local = nr / nprocs
            remainder = mod(nr, nprocs)

            if (rank < remainder) nr_local = nr_local + 1 

            if (rank < remainder) then
                nr_start = rank * nr_local
            else
                nr_start = rank * nr_local + remainder
            end if

            left_rank = rank - 1
            right_rank = rank + 1
            if (rank == 0) left_rank = MPI_PROC_NULL
            if (rank == nprocs - 1) right_rank = MPI_PROC_NULL

            ! Allocate (local + 2 halo columns)
            allocate(phi(0:nz, 0:nr_local+1))
            phi = 0.0_dp

            ! Apply boundary conditions
            do k_local = 1, nr_local
                k_global = nr_start + k_local
                r_val = real(k_global, dp) * p%h

                ! Inner conductor (r <= R_inner)
                if (k_global <= nr_in) then 
                    do i = 0, nz_in
                        phi(i, k_local) = p%phi_inner
                    end do
                end if

                ! Inlet BC (z=0) between inner/outer radius
                if (k_global > nr_in) then
                    phi(0, k_local) = coax_analytic(r_val, p)
                end if
            end do

            ! Axis BC (r=0) on rank 0 only 
            if (rank ==0) then 
                phi(0:nz_in, 0) = p%phi_inner 
            end if 

        end subroutine init_grid

        subroutine probe_point(r_target, z_target, name, p)
            real(kind=dp), intent(in) :: r_target, z_target
            character(len=*), intent(in) :: name
            type(params_t), intent(in) :: p
            integer :: iz, ir, local_r
            real(kind=dp) :: local_phi, global_phi

            iz = nint(z_target / p%h)   ! z-index
            ir = nint(r_target / p%h)   ! r-index

            ! Default invalid value
            local_phi = -1.0_dp      

            ! Check axis (Only Rank 0 has this)
            if (ir == 0) then
                if (rank == 0) local_phi = phi(iz, 0)

            ! Check bulk (Find the rank who owns this column)
            else if (ir > nr_start .and. ir <= nr_start + nr_local) then
                local_r = ir - nr_start
                local_phi = phi(iz, local_r)
            end if
            
            ! REDUCE: Send the highest value found to rank 0 
            ! Ensures order of results
            call MPI_Reduce(local_phi, global_phi, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, comm, ierr)
            
            ! Only rank 0 prints
            if (rank == 0) then
                print '(3A, F10.4, A)', 'Point ', name, ':', global_phi, ' V'
            end if
        end subroutine probe_point      
        
        subroutine save_parallel_dat(p)
            type(params_t), intent(in) :: p
            character(len=64) :: filename
            integer :: unit_num, i, k, k_global
        
            ! Unique filename for each rank
            write(filename, '(A,I0,A)') 'phi_rank_', rank, '.dat'
        
            open(newunit=unit_num, file=filename, status='replace', action='write')
        
            ! Axis (Rank 0 only)
            if (rank == 0) then
                do i = 0, nz
                    ! Z, R, Phi (R=0 here)
                    write(unit_num, *) real(i,dp)*p%h, 0.0_dp, phi(i,0)
                end do
            end if
        
            ! Bulk
            do k = 1, nr_local
                k_global = nr_start + k
                do i = 0, nz
                    write(unit_num, *) real(i,dp)*p%h, real(k_global,dp)*p%h, phi(i,k)
                end do
            end do
        
            close(unit_num)
        
            if (rank == 0) print *, 'Parallel .dat write complete.'
        end subroutine save_parallel_dat

        subroutine cleanup()
            ! Deallocate grid
            if (allocated(phi)) deallocate(phi)
        end subroutine cleanup

end module grid_data

module sor_solver
    use precision
    use sim_params
    use grid_data
    use potentials
    use mpi
    implicit none

contains

    subroutine solve(p)
        type(params_t), intent(in) :: p
        
        integer :: iter, phase, i, k, k_global
        real(kind=dp) :: phi_old, phi_new
        real(kind=dp) :: local_diff, global_diff, diff
        real(kind=dp) :: inv_r_factor
        
        integer :: requests(4)
        integer :: statuses(MPI_STATUS_SIZE, 4)
        integer :: ierr
        
        if (rank == 0) then 
            print *, 'Starting SOR...'
	    print *, repeat('-',30)
        end if 

        do iter = 1, p%max_iter
            local_diff = 0.0_dp

            do phase = 0, 1
                
                ! Start non-blocking communications
                ! Post recvs
                call MPI_Irecv(phi(0, 0), nz+1, MPI_DOUBLE_PRECISION, left_rank, 1, &
                               comm, requests(1), ierr)
                call MPI_Irecv(phi(0, nr_local+1), nz+1, MPI_DOUBLE_PRECISION, right_rank, 2, &
                               comm, requests(2), ierr)
                
                ! Post sends
                call MPI_Isend(phi(0, 1), nz+1, MPI_DOUBLE_PRECISION, left_rank, 2, &
                               comm, requests(3), ierr)
                call MPI_Isend(phi(0, nr_local), nz+1, MPI_DOUBLE_PRECISION, right_rank, 1, &
                               comm, requests(4), ierr)
            
                
                ! Axis update [r=0] (rank 0 only - strictly local)
                if (rank == 0) then
                    k = 0
                    !$OMP PARALLEL DO &
                    !$OMP PRIVATE(i, phi_old, phi_new, diff) &
                    !$OMP REDUCTION(max: local_diff)

                    ! Update along axis
                    do i = nz_in + 1, nz - 1
                        if (mod(i, 2) /= phase) cycle
                        
                        ! SOR update at axis
                        phi_new = update_axis(phi(i,0), phi(i,1), phi(i+1,0), phi(i-1,0), p%omega)
                        
                        ! Store new value and track max difference
                        phi_old = phi(i,0)
                        phi(i, 0) = phi_new
                        diff = abs(phi_new - phi_old)
                        if (diff > local_diff) local_diff = diff
                    end do
                    !$OMP END PARALLEL DO
                end if

                ! Bulk Update
                if (nr_local >= 2) then
                    !$OMP PARALLEL DO & 
                    !$OMP PRIVATE(k, k_global, inv_r_factor, i, phi_old, phi_new, diff) &
                    !$OMP REDUCTION(max: local_diff) &
                    !$OMP SCHEDULE(DYNAMIC)

                    ! Update internal columns
                    do k = 2, nr_local - 1
                        k_global = nr_start + k
                        
                        ! Optimisation: Compute division once per column, pass as multiplier
                        inv_r_factor = 1.0_dp / (8.0_dp * real(k_global, dp))
                        
                        ! Update along z
                        do i = 1, nz - 1
                            ! Skip inner conductor region
                            if (i <= nz_in .and. k_global <= nr_in) cycle

                            ! Red-Black ordering
                            if (mod(i + k_global, 2) /= phase) cycle

                            ! SOR update at bulk
                            phi_new = update_bulk(phi(i,k), phi(i,k-1), phi(i,k+1), &
                                                  phi(i+1,k), phi(i-1,k), &
                                                  inv_r_factor, p%omega)
                            
                            ! Store new value and track max difference
                            phi_old = phi(i,k)
                            phi(i, k) = phi_new
                            diff = abs(phi_new - phi_old)
                            if (diff > local_diff) local_diff = diff
                        end do
                    end do
                    !$OMP END PARALLEL DO
                end if

                ! Wait for comms to complete
                call MPI_Waitall(4, requests, statuses, ierr)


                ! Halo Update
                !$OMP PARALLEL DO & 
                !$OMP PRIVATE(k, k_global, inv_r_factor, i, phi_old, phi_new, diff) &
                !$OMP REDUCTION(max: local_diff)
                do k = 1, nr_local, max(1, nr_local - 1) 
                    if (k > 1 .and. k < nr_local) cycle 

                    ! Determine global r-index
                    k_global = nr_start + k
                    if (k_global >= nr) cycle

                    ! Optimisation: Compute division once per column, pass as multiplier
                    inv_r_factor = 1.0_dp / (8.0_dp * real(k_global, dp))
                    
                    ! Update along z
                    do i = 1, nz - 1
                        if (i <= nz_in .and. k_global <= nr_in) cycle
                        if (mod(i + k_global, 2) /= phase) cycle

                        ! SOR update at bulk
                        phi_new = update_bulk(phi(i,k), phi(i,k-1), phi(i,k+1), &
                                              phi(i+1,k), phi(i-1,k), &
                                              inv_r_factor, p%omega)
                        
                        ! Store new value and track max difference
                        phi_old = phi(i,k)
                        phi(i, k) = phi_new
                        diff = abs(phi_new - phi_old)
                        if (diff > local_diff) local_diff = diff
                    end do
                end do
                !$OMP END PARALLEL DO

            end do

            ! Global convergence
            call MPI_Allreduce(local_diff, global_diff, 1, MPI_DOUBLE_PRECISION, &
                               MPI_MAX, comm, ierr)
	    if (rank == 0) then
                if (iter == 1 .or. mod(iter, 500) == 0) then
                    ! Format: 'CONV' tag, Iteration, Error
                    print '(A, I8, E14.6)', 'CONV ', iter, global_diff
                end if
            end if		

            ! Exit if converged
            if (global_diff < p%tol) exit
        end do
        
        ! Final output
        if (rank == 0) then
            print *, repeat('-',30)
            print *, 'Converged in', iter, 'iterations. Final error:', global_diff
        end if
    end subroutine solve

end module sor_solver

program main
    use precision
    use sim_params
    use grid_data
    use potentials
    use sor_solver
    use mpi
    use omp_lib
    implicit none

    type(params_t) :: params
    real(kind=dp) :: t_start, t_end

    call MPI_Init(ierr)
    
    ! Set up params
    params = default_params()
    
    ! User overrides (for example)
    params%h        = 0.002_dp
    params%omega    = 1.98_dp
    params%tol      = 1.0e-4_dp
    params%max_iter = 1000000
   
    ! Initialise the grid
    call init_grid(params, MPI_COMM_WORLD)

    ! Print configuration (Only Rank 0 prints)
    if (rank == 0) then
        print *, repeat('=', 30)
        print *, 'SIMULATION PARAMETERS'
        print *, repeat('=', 30)
        print *, 'Geometry:'
        print *, '  L_outer:  ', params%L_outer
        print *, '  L_inner:  ', params%L_inner
        print *, '  R_outer:  ', params%R_outer
        print *, '  R_inner:  ', params%R_inner
        print *, '  phi_inner:', params%phi_inner, ' V'
        print *
        print *, 'Solver Settings:'
        print *, '  h:        ', params%h
        print *, '  omega:    ', params%omega
        print *, '  tol:      ', params%tol 
        print *, '  max_iter: ', params%max_iter
        
        print *, repeat('-', 30)
        print *, 'Parallel Config:'
        print *, '  MPI Ranks:     ', nprocs
        print *, '  OpenMP Threads:', omp_get_max_threads()
        print *, '  Global Grid:   ', nz, 'x', nr
        print *, '  Total Points:  ', nz * nr
        print *, repeat('=', 30)
    endif

    ! Run solver
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    t_start = MPI_Wtime()
    call solve(params)
    t_end = MPI_Wtime()

    if (rank == 0) print *, 'Solver Time:', t_end - t_start, 's'

    ! Check results
    if (rank == 0) print *, '--- Final Results ---'
    call probe_point(0.0_dp,  7.5_dp,  'A', params)
    call probe_point(6.0_dp,  5.0_dp,  'B', params)
    call probe_point(5.0_dp, 12.5_dp, 'C', params)
    
    ! Save the field for plotting
    call save_parallel_dat(params)

    ! Cleanup and close MPI
    call cleanup()
    call MPI_Finalize(ierr)

end program main
