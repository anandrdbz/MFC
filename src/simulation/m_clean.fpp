!>
!! @file m_clean.fpp

module m_clean

    ! Dependencies =============================================================
    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    use m_start_up             !< Reading and checking procedures for the input
    !< and the initial condition and grid data files

    use m_variables_conversion !< State variables type conversion procedures

    use m_weno                 !< Weighted and essentially non-oscillatory (WENO)
                               !! schemes for spatial reconstruction of variables

    use m_riemann_solvers      !< Exact and approximate Riemann problem solvers

    use m_cbc                  !< Characteristic boundary conditions (CBC)

    use m_monopole             !< Monopole calculations

    use m_rhs                  !< Right-hand-side (RHS) evaluation procedures

    use m_data_output          !< Run-time info & solution data output procedures

    use m_time_steppers        !< Time-stepping algorithms

    use m_qbmm                 !< Quadrature MOM

    use m_derived_variables     !< Procedures used to compute quantites derived
                                !! from the conservative and primitive variables

    use m_hypoelastic

    use m_viscous

    use m_bubbles

    use ieee_arithmetic

#ifdef _OPENACC
    use openacc
#endif

    use m_nvtx

    use m_helper
    ! ==========================================================================

    implicit none

    private; public :: s_initialize_modules, s_initialize_gpu_vars, s_initialize_mpi_domain, s_finalize_modules, &
                        s_perform_time_step, s_save_data

contains

    subroutine s_perform_time_step(t_step, time_avg, time_final, io_time_avg, io_time_final, proc_time, io_proc_time, file_exists, start, finish, nt)
        integer, intent(INOUT) :: t_step
        real(kind(0d0)), intent(INOUT) :: time_avg, time_final
        real(kind(0d0)), intent(INOUT) :: io_time_avg, io_time_final
        real(kind(0d0)),  dimension(:), intent(INOUT) :: proc_time
        real(kind(0d0)),  dimension(:), intent(INOUT) :: io_proc_time
        logical, intent(INOUT) :: file_exists
        real(kind(0d0)), intent(INOUT) :: start, finish
        integer, intent(INOUT) :: nt

        integer :: i, j, k, l

        if (proc_rank == 0) then
            print '(" ["I3"%]  Time step "I8" of "I0" @ t_step = "I0"")',                             &
                  int(ceiling(100d0*(real(t_step - t_step_start)/(t_step_stop - t_step_start + 1)))), &
                  t_step      - t_step_start + 1,                                                     &
                  t_step_stop - t_step_start + 1,                                                     &
                  t_step
        end if
        mytime = mytime + dt

        if (probe_wrt) then
            do i = 1, sys_size
                !$acc update host(q_cons_ts(1)%vf(i)%sf)
            end do
        end if

        call s_compute_derived_variables(t_step)

#ifdef DEBUG
        print *, 'Computed derived vars'
#endif

        ! Total-variation-diminishing (TVD) Runge-Kutta (RK) time-steppers
        if (time_stepper == 1) then
            call s_1st_order_tvd_rk(t_step, time_avg)
        elseif (time_stepper == 2) then
            call s_2nd_order_tvd_rk(t_step, time_avg)
        elseif (time_stepper == 3) then
            call s_3rd_order_tvd_rk(t_step, time_avg)
        end if

        ! Time-stepping loop controls

        if (t_step == t_step_stop) then

            call s_mpi_barrier()

            if (num_procs > 1) then
                call mpi_bcast_time_step_values(proc_time, time_avg)

                call mpi_bcast_time_step_values(io_proc_time, io_time_avg)
            end if

            if (proc_rank == 0) then
                time_final = 0d0
                io_time_final = 0d0
                if (num_procs == 1) then
                    time_final = time_avg
                    io_time_final = io_time_avg
                    print *, "Final Time", time_final
                else
                    time_final = maxval(proc_time)
                    io_time_final = maxval(io_proc_time)
                    print *, "Final Time", time_final
                end if
                inquire (FILE='time_data.dat', EXIST=file_exists)
                if (file_exists) then
                    open (1, file='time_data.dat', position='append', status='old')
                    write (1, *) num_procs, time_final
                    close (1)
                else
                    open (1, file='time_data.dat', status='new')
                    write (1, *) num_procs, time_final
                    close (1)
                end if

                inquire (FILE='io_time_data.dat', EXIST=file_exists)
                if (file_exists) then
                    open (1, file='io_time_data.dat', position='append', status='old')
                    write (1, *) num_procs, io_time_final
                    close (1)
                else
                    open (1, file='io_time_data.dat', status='new')
                    write (1, *) num_procs, io_time_final
                    close (1)
                end if

            end if

            
        else
            if ((mytime + dt) >= finaltime) dt = finaltime - mytime
            t_step = t_step + 1
        end if
    end subroutine s_perform_time_step

    subroutine s_save_data(t_step, start, finish, io_time_avg, nt)
        real(kind(0d0)), intent(INOUT) ::  start, finish, io_time_avg
        integer, intent(INOUT) :: t_step, nt
        integer :: i, j, k, l

        if (mod(t_step - t_step_start, t_step_save) == 0 .or. t_step == t_step_stop) then

            call cpu_time(start)
            !  call nvtxStartRange("I/O")
            do i = 1, sys_size
                !$acc update host(q_cons_ts(1)%vf(i)%sf)
                do l = 0, p
                    do k = 0, n
                        do j = 0, m
                            if(ieee_is_nan(q_cons_ts(1)%vf(i)%sf(j, k, l))) then
                                print *, "NaN(s) in timestep output.", j, k, l, i,  proc_rank, t_step, m, n, p                                
                                error stop "NaN(s) in timestep output."
                            end if
                        end do
                    end do
                end do
            end do

            if(qbmm .and. .not. polytropic) then
                !$acc update host(pb_ts(1)%sf)
                !$acc update host(mv_ts(1)%sf)
            end if

            call s_write_data_files(q_cons_ts(1)%vf, q_prim_vf, t_step)
            !  call nvtxEndRange
            call cpu_time(finish)
            nt = int((t_step - t_step_start)/(t_step_save))
            if (nt == 1) then
                io_time_avg = abs(finish - start)
            else
                io_time_avg = (abs(finish - start) + io_time_avg*(nt - 1))/nt
            end if
        end if
    end subroutine s_save_data    

    subroutine s_initialize_modules()
        call s_initialize_global_parameters_module()
        !Quadrature weights and nodes for polydisperse simulations
        if(bubbles .and. nb > 1 .and. R0_type == 1) then
            call s_simpson
        end if
        !Initialize variables for non-polytropic (Preston) model
        if(bubbles .and. .not. polytropic) then
            call s_initialize_nonpoly()
        end if
        !Initialize pb based on surface tension for qbmm (polytropic)
        if(qbmm .and. polytropic .and. Web /= dflt_real) then
            pb0 = pref + 2d0 * fluid_pp(1)%ss / (R0*R0ref)
            pb0 = pb0 / pref
            pref = 1d0                              
        end if

#if defined(_OPENACC) && defined(MFC_MEMORY_DUMP)
        call acc_present_dump()
#endif

        call s_initialize_mpi_proxy_module()
        call s_initialize_variables_conversion_module()
        if (grid_geometry == 3) call s_initialize_fftw_module()
        call s_initialize_start_up_module()
        call s_initialize_riemann_solvers_module()

        if(bubbles) call s_initialize_bubbles_module()

        if (qbmm) call s_initialize_qbmm_module()

#if defined(_OPENACC) && defined(MFC_MEMORY_DUMP)
        call acc_present_dump()
#endif

        if (monopole) then
            call s_initialize_monopole_module()
        end if
        if (any(Re_size > 0)) then
            call s_initialize_viscous_module()
        end if
        call s_initialize_rhs_module()

#if defined(_OPENACC) && defined(MFC_MEMORY_DUMP)
        call acc_present_dump()
#endif

        if (hypoelasticity) call s_initialize_hypoelastic_module()
        call s_initialize_data_output_module()
        call s_initialize_derived_variables_module()
        call s_initialize_time_steppers_module()

#if defined(_OPENACC) && defined(MFC_MEMORY_DUMP)
        call acc_present_dump()
#endif

        ! Associate pointers for serial or parallel I/O
        if (parallel_io .neqv. .true.) then
            s_read_data_files => s_read_serial_data_files
            s_write_data_files => s_write_serial_data_files
        else
            s_read_data_files => s_read_parallel_data_files
            s_write_data_files => s_write_parallel_data_files
        end if

        ! Reading in the user provided initial condition and grid data
        call s_read_data_files(q_cons_ts(1)%vf)
        if (model_eqns == 3) call s_initialize_internal_energy_equations(q_cons_ts(1)%vf)

        ! Populating the buffers of the grid variables using the boundary conditions
        call s_populate_grid_variables_buffers()

        ! Computation of parameters, allocation of memory, association of pointers,
        ! and/or execution of any other tasks that are needed to properly configure
        ! the modules. The preparations below DO DEPEND on the grid being complete.
        call s_initialize_weno_module()

#if defined(_OPENACC) && defined(MFC_MEMORY_DUMP)
        print *, "[MEM-INST] After: s_initialize_weno_module"
        call acc_present_dump()
#endif

        call s_initialize_cbc_module()

        call s_initialize_derived_variables()

    end subroutine s_initialize_modules

    subroutine s_initialize_mpi_domain()
        integer :: ierr
#ifdef _OPENACC
        real(kind(0d0)) :: starttime, endtime
        integer :: num_devices, local_size, num_nodes, ppn, my_device_num
        integer :: dev, devNum, local_rank
#ifdef MFC_MPI
        integer :: local_comm
#endif
        integer(acc_device_kind) :: devtype
#endif

    ! Initializing MPI execution environment

        call s_mpi_initialize()

    ! Bind GPUs if OpenACC is enabled
#ifdef _OPENACC
#ifndef MFC_MPI
        local_size = 1
        local_rank = 0
#else
        call MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, &
                                 MPI_INFO_NULL, local_comm, ierr)
        call MPI_Comm_size(local_comm, local_size, ierr)
        call MPI_Comm_rank(local_comm, local_rank, ierr)
#endif

        devtype = acc_get_device_type()
        devNum  = acc_get_num_devices(devtype)
        dev     = mod(local_rank, devNum)

        call acc_set_device_num(dev, devtype)
#endif

        ! The rank 0 processor assigns default values to the user inputs prior to
        ! reading them in from the input file. Next, the user inputs are read and
        ! their consistency is checked. The identification of any inconsistencies
        ! will result in the termination of the simulation.
        if (proc_rank == 0) then
            call s_assign_default_values_to_user_inputs()
            call s_read_input_file()
            call s_check_input_file()
            print '(" Simulating a "I0"x"I0"x"I0" case on "I0" rank(s)")', m, n, p, num_procs
        end if

        ! Broadcasting the user inputs to all of the processors and performing the
        ! parallel computational domain decomposition. Neither procedure has to be
        ! carried out if the simulation is in fact not truly executed in parallel.

        call s_mpi_bcast_user_inputs()
        call s_initialize_parallel_io()
        call s_mpi_decompose_computational_domain()

    end subroutine s_initialize_mpi_domain

    subroutine s_initialize_gpu_vars()
        integer :: i
        !Update GPU DATA
        !$acc update device(dt, dx, dy, dz, x_cc, y_cc, z_cc, x_cb, y_cb, z_cb)
        !$acc update device(sys_size, buff_size)
        !$acc update device(m, n, p)
        !$acc update device(momxb, momxe, bubxb, bubxe, advxb, advxe, contxb, contxe, strxb, strxe)
        do i = 1, sys_size
        !$acc update device(q_cons_ts(1)%vf(i)%sf)
        end do
        if(qbmm .and. .not. polytropic) then
        !$acc update device(pb_ts(1)%sf, mv_ts(1)%sf)
        end if
        !$acc update device(dt, sys_size, pref, rhoref, gamma_idx, pi_inf_idx, E_idx, alf_idx, stress_idx, mpp_lim, bubbles, hypoelasticity, alt_soundspeed, avg_state, num_fluids, model_eqns, num_dims, mixture_err, nb, weight, grid_geometry, cyl_coord, mapped_weno, mp_weno, weno_eps)
        !$acc update device(nb, R0ref, Ca, Web, Re_inv, weight, R0, V0, bubbles, polytropic, polydisperse, qbmm, R0_type, ptil, bubble_model, thermal, poly_sigma)
        !$acc update device(R_n, R_v, phi_vn, phi_nv, Pe_c, Tw, pv, M_n, M_v, k_n, k_v, pb0, mass_n0, mass_v0, Pe_T, Re_trans_T, Re_trans_c, Im_trans_T, Im_trans_c, omegaN , mul0, ss, gamma_v, mu_v, gamma_m, gamma_n, mu_n, gam)
        !$acc update device(monopole, num_mono)
    end subroutine s_initialize_gpu_vars


    subroutine s_finalize_modules()
        ! Disassociate pointers for serial and parallel I/O
        s_read_data_files => null()
        s_write_data_files => null()

        call s_finalize_time_steppers_module()
        call s_finalize_derived_variables_module()
        call s_finalize_data_output_module()
        call s_finalize_rhs_module()
        call s_finalize_cbc_module()
        call s_finalize_riemann_solvers_module()
        call s_finalize_weno_module()
        call s_finalize_start_up_module()
        call s_finalize_variables_conversion_module()
        if (grid_geometry == 3) call s_finalize_fftw_module
        call s_finalize_mpi_proxy_module()
        call s_finalize_global_parameters_module()

        if (any(Re_size > 0)) then
            call s_finalize_viscous_module()
        end if

        ! Terminating MPI execution environment
        call s_mpi_finalize()
    end subroutine s_finalize_modules

end module m_clean