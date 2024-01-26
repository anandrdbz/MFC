!>
!! @file m_assign_variables.f90
!! @brief Contains module m_assign_variables
module m_assign_variables

    ! Dependencies =============================================================
    use m_derived_types         ! Definitions of the derived types

    use m_global_parameters     ! Global parameters for the code

    use m_variables_conversion  ! Subroutines to change the state variables from
    ! one form to another
    ! ==========================================================================

    implicit none

    public :: s_perturb_primitive

    type(scalar_field) :: alf_sum

    procedure(s_assign_patch_xxxxx_primitive_variables), &
    pointer :: s_assign_patch_primitive_variables => null() !<
    !! Depending on the multicomponent flow model, this variable is a pointer to
    !! either the subroutine s_assign_patch_mixture_primitive_variables, or the
    !! subroutine s_assign_patch_species_primitive_variables

    !> Abstract interface to the two subroutines that assign the patch primitive
    !! variables, either mixture or species, depending on the subroutine, to a
    !! particular cell in the computational domain
    abstract interface

        !> Skeleton of s_assign_patch_mixture_primitive_variables
        !!      and s_assign_patch_species_primitive_variables
        !! @param patch_id is the patch identifier
        !! @param j (x) cell index in which the mixture or species primitive variables from the indicated patch areassigned
        !! @param k (y,th) cell index in which the mixture or species primitive variables from the indicated patch areassigned
        !! @param l (z) cell index in which the mixture or species primitive variables from the indicated patch areassigned
        subroutine s_assign_patch_xxxxx_primitive_variables(patch_id, j, k, l, &
                                                eta, q_prim_vf, patch_id_fp)

            import :: scalar_field, sys_size, n, m, p

            integer, intent(IN) :: patch_id
            integer, intent(IN) :: j, k, l
            integer, intent(INOUT), dimension(0:m, 0:n, 0:p) :: patch_id_fp
            type(scalar_field), dimension(1:sys_size) :: q_prim_vf
            real(kind(0d0)) :: eta !<

        end subroutine s_assign_patch_xxxxx_primitive_variables

    end interface

    private; public :: s_initialize_assign_variables_module, &
        s_assign_patch_primitive_variables, &
        s_assign_patch_mixture_primitive_variables, &
        s_assign_patch_species_primitive_variables, &
        s_finialize_assign_variables_module

contains

    subroutine s_initialize_assign_variables_module()

        allocate (alf_sum%sf(0:m, 0:n, 0:p))

        ! Depending on multicomponent flow model, the appropriate procedure
        ! for assignment of the patch mixture or species primitive variables
        ! to a cell in the domain is targeted by the procedure pointer

        if (model_eqns == 1) then        ! Gamma/pi_inf model
            s_assign_patch_primitive_variables => &
                s_assign_patch_mixture_primitive_variables
        else ! Volume fraction model
            s_assign_patch_primitive_variables => &
                s_assign_patch_species_primitive_variables
        end if
    
    end subroutine s_initialize_assign_variables_module

    !>  This subroutine assigns the mixture primitive variables
        !!              of the patch designated by the patch_id, to the cell that
        !!              is designated by the indexes (j,k,l). In addition, the
        !!              variable bookkeeping the patch identities in the entire
        !!              domain is updated with the new assignment. Note that if
        !!              the smoothing of the patch's boundaries is employed, the
        !!              ensuing primitive variables in the cell will be a type of
        !!              combination of the current patch's primitive variables
        !!              with those of the smoothing patch. The specific details
        !!              of the combination may be found in Shyue's work (1998).
        !! @param patch_id the patch identifier
        !! @param j  the x-dir node index
        !! @param k  the y-dir node index
        !! @param l  the z-dir node index
    subroutine s_assign_patch_mixture_primitive_variables(patch_id, j, k, l, &
                                         eta, q_prim_vf, patch_id_fp)

        !$acc routine seq
        integer, intent(IN) :: patch_id
        integer, intent(INOUT), dimension(0:m, 0:n, 0:p) :: patch_id_fp
        type(scalar_field), dimension(1:sys_size) :: q_prim_vf
        real(kind(0d0)) :: eta !<

        integer, intent(IN) :: j, k, l

        real(kind(0d0)) :: rho    !< density
        real(kind(0d0)), dimension(int(E_idx - mom_idx%beg)) :: vel    !< velocity
        real(kind(0d0)) :: pres   !< pressure
        real(kind(0d0)) :: gamma  !< specific heat ratio function
        real(kind(0d0)) :: x_centroid, y_centroid
        real(kind(0d0)) :: epsilon, beta

        integer :: smooth_patch_id
        integer :: i !< generic loop operator

        ! Assigning the mixture primitive variables of a uniform state patch

        ! Transferring the identity of the smoothing patch
        smooth_patch_id = patch_icpp(patch_id)%smooth_patch_id

        ! Density
        q_prim_vf(1)%sf(j, k, l) = &
            eta*patch_icpp(patch_id)%rho &
            + (1d0 - eta)*patch_icpp(smooth_patch_id)%rho

        ! Velocity
        do i = 1, E_idx - mom_idx%beg
            q_prim_vf(i + 1)%sf(j, k, l) = &
                1d0/q_prim_vf(1)%sf(j, k, l)* &
                (eta*patch_icpp(patch_id)%rho &
                    *patch_icpp(patch_id)%vel(i) &
                    + (1d0 - eta)*patch_icpp(smooth_patch_id)%rho &
                    *patch_icpp(smooth_patch_id)%vel(i))
        end do

        ! Specific heat ratio function
        q_prim_vf(gamma_idx)%sf(j, k, l) = &
            eta*patch_icpp(patch_id)%gamma &
            + (1d0 - eta)*patch_icpp(smooth_patch_id)%gamma

        ! Pressure
        q_prim_vf(E_idx)%sf(j, k, l) = &
            1d0/q_prim_vf(gamma_idx)%sf(j, k, l)* &
            (eta*patch_icpp(patch_id)%gamma &
                *patch_icpp(patch_id)%pres &
                + (1d0 - eta)*patch_icpp(smooth_patch_id)%gamma &
                *patch_icpp(smooth_patch_id)%pres)

        ! Liquid stiffness function
        q_prim_vf(pi_inf_idx)%sf(j, k, l) = &
            eta*patch_icpp(patch_id)%pi_inf &
            + (1d0 - eta)*patch_icpp(smooth_patch_id)%pi_inf

        ! Updating the patch identities bookkeeping variable
        if (1d0 - eta < 1d-16) patch_id_fp(j, k, l) = patch_id

    end subroutine s_assign_patch_mixture_primitive_variables ! ------------

    function f_bump(j, k)
        integer, intent(IN) :: j, k
        real(kind(0d0)) :: f_bump
        real(kind(0d0)) :: r 

        r = sqrt((x_cc(j) - 0.5)**2d0 + (y_cc(k) - 0.5)**2d0)
        if(r >= 0.3d0) then
            f_bump = 0d0
        else
            f_bump = exp(-0.01d0 / (0.3d0**2d0 - r**2d0))
        end if


    end function f_bump

    !Stable perturbation in pressure (Ando)
    subroutine s_perturb_primitive(j, k, l, q_prim_vf)

        type(scalar_field), dimension(1:sys_size), intent(INOUT) :: q_prim_vf
        integer, intent(IN) :: j, k, l

        integer :: i
        real(kind(0d0)) :: pres_mag , loc, n_tait, B_tait, p0
        real(kind(0d0)) :: R3bar, n0, ratio, nH, vfH, velH, rhoH, deno
        real(kind(0d0)) :: x, y
        real(kind(0d0)) :: pi = 3.141563

        p0 = 101325
        pres_mag = 1D-1
        !loc = x_cc(177)
        n_tait = fluid_pp(1)%gamma
        B_tait = fluid_pp(1)%pi_inf

        n_tait = 1.d0/n_tait + 1.d0 
        B_tait = B_tait * (n_tait - 1d0) / n_tait

        ! if(x_cc(j) >= 0.25d0 .and. x_cc(j) <= 0.75d0) then
        !     q_prim_vf(momxb)%sf(j, k, l) = -1d0*sin(2*3.141563*(x_cc(j) - 0.5d0)/1.d0)
        ! else if(x_cc(j) < 0.25) then
        !     q_prim_vf(momxb)%sf(j, k, l) = -1d0*sin(2*3.141563*(0.25d0 - 0.5d0)/1.d0)
        ! else if(x_cc(j) > 0.75) then
        !     q_prim_vf(momxb)%sf(j, k, l) = -1d0*sin(2*3.141563*(0.75d0 - 0.5d0)/1.d0)
        ! end if

        ! if(y_cc(k) >= 0.25d0 .and. y_cc(k) <= 0.75d0) then        
        !     q_prim_vf(momxb + 1)%sf(j, k, l) = -1d0*sin(2*3.141563*(y_cc(k) - 0.5d0)/1.d0)
        ! else if(y_cc(k) < 0.25) then
        !     q_prim_vf(momxb + 1)%sf(j, k, l) = -1d0*sin(2*3.141563*(0.25d0 - 0.5d0)/1.d0)
        ! else if(y_cc(k) > 0.75) then
        !     q_prim_vf(momxb + 1)%sf(j, k, l) = -1d0*sin(2*3.141563*(0.75d0 - 0.5d0)/1.d0)
        ! end if



        ! if(j > k) then
        !     q_prim_vf(momxb)%sf(j, k, l) = 1d0*sin(2*3.141563*x_cc(j))*f_bump(j, k) 
        ! else if(k > j) then
        !     q_prim_vf(momxb)%sf(j, k, l) = 1d0*sin(2*3.141563*(x_cc(j) - 0.5))*f_bump(j, k) 
        ! else
        !     q_prim_vf(momxb)%sf(j, k, l) = 0d0
        ! end if

        ! if(j > k) then
        !     q_prim_vf(momxb + 1)%sf(j, k, l) = 1d0*sin(2*3.141563*y_cc(k))*f_bump(j, k) 
        ! else if(k > j) then
        !     q_prim_vf(momxb + 1)%sf(j, k, l) = 1d0*sin(2*3.141563*(y_cc(k) - 0.5))*f_bump(j, k) 
        ! else
        !     q_prim_vf(momxb + 1)%sf(j, k, l) = 0d0
        ! end if

        q_prim_vf(1)%sf(j, k, l) = 1d0
        q_prim_vf(1)%sf(j, k, l) = q_prim_vf(1)%sf(j, k, l) + (0.25d0 / (3d0*0.08d0)) * exp(-0.5d0*( (x_cc(j) - 0.4d0)**2d0 + (y_cc(k) - 0.45d0)**2d0)/ 0.08d0**2d0) 
        q_prim_vf(1)%sf(j, k, l) = q_prim_vf(1)%sf(j, k, l) + (0.25d0 / (3d0*0.08d0)) * exp(-0.5d0*( (x_cc(j) - 0.55d0)**2d0 + (y_cc(k) - 0.55d0)**2d0)/ 0.08d0**2d0) 
        q_prim_vf(1)%sf(j, k, l) = q_prim_vf(1)%sf(j, k, l) + (0.25d0 / (3d0*0.08d0)) * exp(-0.5d0*( (x_cc(j) - 0.6d0)**2d0 + (y_cc(k) - 0.42d0)**2d0)/ 0.08d0**2d0) 
        q_prim_vf(E_idx)%sf(j, k, l) = 0.2 * (q_prim_vf(1)%sf(j, k, l))**1.4d0

        ! x = x_cc(j) - 0.5
        ! y = y_cc(k) - 0.5
        ! if(x**2d0 + y**2d0 <= 0.125d0**2d0) then
        !     q_prim_vf(momxb)%sf(j, k, l) = 1d0 * x / sqrt(x**2d0 + y**2) !* tanh((0.125d0 - sqrt(x**2 + y**2)))
        !     q_prim_vf(momxb + 1)%sf(j, k, l) = 1d0 * y / sqrt(x**2d0 + y**2) !* tanh((0.125d0 - sqrt(x**2 + y**2)))
        ! end if

        ! if(x**2d0 + y**2d0 <= 0.125d0**2d0) then
        !     q_prim_vf(momxb)%sf(j, k, l) = 1d0 * x / sqrt(x**2d0 + y**2) !* tanh((0.125d0 - sqrt(x**2 + y**2)))
        !     q_prim_vf(momxb + 1)%sf(j, k, l) = 1d0 * y / sqrt(x*2d0 + y**2) !* tanh((0.125d0 - sqrt(x**2 + y**2)))
        ! end if

        ! if( (x**2d0 + y**2d0 >0.125d0**2d0) .and. (x**2d0 + y**2d0 <= 0.25d0**2d0) ) then
        !     q_prim_vf(momxb)%sf(j, k, l) = -1d0 * x / sqrt(x**2d0 + y**2) !* tanh((0.125d0 - sqrt(x**2 + y**2)))
        !     q_prim_vf(momxb + 1)%sf(j, k, l) = -1d0 * y / sqrt(x*2d0 + y**2) !* tanh((0.125d0 - sqrt(x**2 + y**2)))
        ! end if

#if 0
! C         if(j < 177) then
! C             q_prim_vf(E_idx)%sf(j, k, l) = 0.5 * q_prim_vf(E_idx)%sf(j, k, l) 
! C         end if


! C         if(qbmm) then
! C             do i = 1, nb
! C                 q_prim_vf(bubxb + 1 + (i-1)*nmom)%sf(j, k, l) = q_prim_vf(bubxb + 1 + (i-1)*nmom)%sf(j, k, l) * ((p0 - fluid_pp(1)%pv) / (q_prim_vf(E_idx)%sf(j, k, l) * p0 - fluid_pp(1)%pv)) ** (1 / 3d0)
! C             end do
! C         end if


! C         R3bar = 0d0

! C         if(qbmm) then
! C             do i = 1, nb
! C                 R3bar = R3bar + weight(i) * 0.5d0 * (q_prim_vf(bubxb + 1 + (i-1)*nmom)%sf(j, k, l) ) ** 3d0
! C                 R3bar = R3bar + weight(i) * 0.5d0 * (q_prim_vf(bubxb + 1 + (i-1)*nmom)%sf(j, k, l) ) ** 3d0
! C             end do
! C         else
! C             do i = 1, nb
! C                 if(polytropic) then
! C                     R3bar = R3bar + weight(i) * (q_prim_vf(bubxb + (i - 1) * 2)%sf(j, k, l)) ** 3d0
! C                 else
! C                     R3bar = R3bar + weight(i) * (q_prim_vf(bubxb + (i - 1) * 4)%sf(j, k, l)) ** 3d0
! C                 end if
! C             end do
! C         end if

! C         n0 = 3d0 * q_prim_vf(alf_idx) % sf(j, k, l) / (4d0 * pi * R3bar)

! C         ratio = ((1d0 + B_tait) / (q_prim_vf(E_idx)%sf(j, k, l) + B_tait)) ** (1D0 / n_tait)

! C         nH = n0 / ( (1d0 - q_prim_vf(alf_idx)%sf(j, k, l)) * ratio + (4d0 * pi / 3d0) * n0 * R3bar )
! C         vfH = (4d0 * pi / 3d0) * nH * R3bar
! C         rhoH = (1d0 - vfH) / ratio
! C         deno = 1d0 - (1d0 - q_prim_vf(alf_idx)%sf(j, k, l)) / rhoH

! C         if(deno == 0d0) then
! C             velH = 0d0
! C         else
! C             velH = (q_prim_vf(E_idx)%sf(j, k, l) - 1d0) / (1d0 - q_prim_vf(alf_idx)%sf(j, k, l)) / deno
! C             velH = dsqrt(velH)
! C             velH = velH * deno
! C         end if

! C         do i = cont_idx%beg, cont_idx%end
! C             q_prim_vf(i)%sf(j, k, l) = rhoH
! C         end do

! C         do i = mom_idx%beg, mom_idx%end
! C             q_prim_vf(i)%sf(j, k, l) = velH
! C         end do

! C         q_prim_vf(alf_idx)%sf(j, k, l) = vfH
#endif
    end subroutine s_perturb_primitive

    !>  This subroutine assigns the species primitive variables. This follows
        !!  s_assign_patch_species_primitive_variables with adaptation for
        !!  ensemble-averaged bubble modeling
        !! @param patch_id the patch identifier
        !! @param j  the x-dir node index
        !! @param k  the y-dir node index
        !! @param l  the z-dir node index
    subroutine s_assign_patch_species_primitive_variables(patch_id, j, k, l, &
                                                eta, q_prim_vf, patch_id_fp)

        !$acc routine seq
        integer, intent(IN) :: patch_id
        integer, intent(INOUT), dimension(0:m, 0:n, 0:p) :: patch_id_fp
        type(scalar_field), dimension(1:sys_size) :: q_prim_vf
        real(kind(0d0)) :: eta !<
        integer, intent(IN) :: j, k, l

        ! Density, the specific heat ratio function and the liquid stiffness
        ! function, respectively, obtained from the combination of primitive
        ! variables of the current and smoothing patches
        real(kind(0d0)) :: rho          !< density
        real(kind(0d0)) :: gamma
        real(kind(0d0)) :: lit_gamma    !< specific heat ratio
        real(kind(0d0)) :: pi_inf       !< stiffness from SEOS
        real(kind(0d0)) :: orig_rho
        real(kind(0d0)) :: orig_gamma
        real(kind(0d0)) :: orig_pi_inf
        real(kind(0d0)) :: muR, muV

        real(kind(0d0)), dimension(int(E_idx - mom_idx%beg)) :: vel    !< velocity
        real(kind(0d0)) :: pres   !< pressure
        real(kind(0d0)) :: x_centroid, y_centroid
        real(kind(0d0)) :: epsilon, beta

        real(kind(0d0)), dimension(sys_size) :: orig_prim_vf !<
            !! Vector to hold original values of cell for smoothing purposes

        integer :: i  !< Generic loop iterator
        integer :: smooth_patch_id

        ! Transferring the identity of the smoothing patch
        smooth_patch_id = patch_icpp(patch_id)%smooth_patch_id
        ! print*, 'smoothing patch id', smooth_patch_id

        ! Transferring original primitive variables
        do i = 1, sys_size
            orig_prim_vf(i) = q_prim_vf(i)%sf(j, k, l)
        end do

        if (mpp_lim .and. bubbles) then
            !adjust volume fractions, according to modeled gas void fraction
            alf_sum%sf = 0d0
            do i = adv_idx%beg, adv_idx%end - 1
                alf_sum%sf = alf_sum%sf + q_prim_vf(i)%sf
            end do

            do i = adv_idx%beg, adv_idx%end - 1
                q_prim_vf(i)%sf = q_prim_vf(i)%sf*(1.d0 - q_prim_vf(alf_idx)%sf) &
                                  /alf_sum%sf
            end do
        end if

        ! Computing Mixture Variables from Original Primitive Variables
        ! call s_convert_species_to_mixture_variables( &
        ! print*, 'first call to convert to mixture'
        ! call s_convert_to_mixture_variables( &
        !     q_prim_vf, j, k, l, &
        !     orig_rho, &
        !     orig_gamma, &
        !     orig_pi_inf)
            ! patch_icpp(patch_id)%rho, &
            ! patch_icpp(patch_id)%gamma, &
            ! patch_icpp(patch_id)%pi_inf)


        ! Computing Mixture Variables of Current Patch =====================

        ! Volume fraction(s)
        do i = adv_idx%beg, adv_idx%end
            q_prim_vf(i)%sf(j, k, l) = patch_icpp(patch_id)%alpha(i - E_idx)
        end do

        call s_convert_to_mixture_variables( &
            q_prim_vf, j, k, l, &
            patch_icpp(patch_id)%rho, &
            patch_icpp(patch_id)%gamma, &
            patch_icpp(patch_id)%pi_inf)
            ! orig_rho, &
            ! orig_gamma, &
            ! orig_pi_inf)

        if (mpp_lim .and. bubbles) then
            !adjust volume fractions, according to modeled gas void fraction
            alf_sum%sf = 0d0
            do i = adv_idx%beg, adv_idx%end - 1
                alf_sum%sf = alf_sum%sf + q_prim_vf(i)%sf
            end do

            do i = adv_idx%beg, adv_idx%end - 1
                q_prim_vf(i)%sf = q_prim_vf(i)%sf*(1.d0 - q_prim_vf(alf_idx)%sf) &
                                  /alf_sum%sf
            end do
        end if

        ! Partial densities
        if (model_eqns /= 4) then
            do i = 1, cont_idx%end
                q_prim_vf(i)%sf(j, k, l) = patch_icpp(patch_id)%alpha_rho(i)
            end do
        end if
  
        ! Density and the specific heat ratio and liquid stiffness functions
        ! call s_convert_species_to_mixture_variables( &
        ! print*, 'second call to convert to mixture'
        call s_convert_to_mixture_variables( &
            q_prim_vf, j, k, l, &
            patch_icpp(patch_id)%rho, &
            patch_icpp(patch_id)%gamma, &
            patch_icpp(patch_id)%pi_inf)

        ! ==================================================================

        ! Computing Mixture Variables of Smoothing Patch ===================

        if (model_eqns /= 4) then
            ! Partial densities
            do i = 1, cont_idx%end
                q_prim_vf(i)%sf(j, k, l) = patch_icpp(smooth_patch_id)%alpha_rho(i)
            end do
        end if

        ! Volume fraction(s)
        do i = adv_idx%beg, adv_idx%end
            q_prim_vf(i)%sf(j, k, l) = patch_icpp(smooth_patch_id)%alpha(i - E_idx)
        end do

        if (mpp_lim .and. bubbles) then
            !adjust volume fractions, according to modeled gas void fraction
            alf_sum%sf = 0d0
            do i = adv_idx%beg, adv_idx%end - 1
                alf_sum%sf = alf_sum%sf + q_prim_vf(i)%sf
            end do

            do i = adv_idx%beg, adv_idx%end - 1
                q_prim_vf(i)%sf = q_prim_vf(i)%sf*(1.d0 - q_prim_vf(alf_idx)%sf) &
                                  /alf_sum%sf
            end do
        end if

        ! Bubbles variables
        if (bubbles) then
            do i = 1, nb
                muR = R0(i)*patch_icpp(smooth_patch_id)%r0 ! = R0(i)
                muV = V0(i)*patch_icpp(smooth_patch_id)%v0 ! = 0
                if (qbmm) then
                    ! Initialize the moment set
                    if (dist_type == 1) then
                        q_prim_vf(bub_idx%fullmom(i, 0, 0))%sf(j, k, l) = 1d0
                        q_prim_vf(bub_idx%fullmom(i, 1, 0))%sf(j, k, l) = muR
                        q_prim_vf(bub_idx%fullmom(i, 0, 1))%sf(j, k, l) = muV
                        q_prim_vf(bub_idx%fullmom(i, 2, 0))%sf(j, k, l) = muR**2 + sigR**2
                        q_prim_vf(bub_idx%fullmom(i, 1, 1))%sf(j, k, l) = muR*muV + rhoRV*sigR*sigV
                        q_prim_vf(bub_idx%fullmom(i, 0, 2))%sf(j, k, l) = muV**2 + sigV**2
                    else if (dist_type == 2) then
                        q_prim_vf(bub_idx%fullmom(i, 0, 0))%sf(j, k, l) = 1d0
                        q_prim_vf(bub_idx%fullmom(i, 1, 0))%sf(j, k, l) = dexp((sigR**2)/2d0)*muR
                        q_prim_vf(bub_idx%fullmom(i, 0, 1))%sf(j, k, l) = muV
                        q_prim_vf(bub_idx%fullmom(i, 2, 0))%sf(j, k, l) = dexp((sigR**2)*2d0)*(muR**2)
                        q_prim_vf(bub_idx%fullmom(i, 1, 1))%sf(j, k, l) = dexp((sigR**2)/2d0)*muR*muV
                        q_prim_vf(bub_idx%fullmom(i, 0, 2))%sf(j, k, l) = muV**2 + sigV**2
                    end if
                else
                    q_prim_vf(bub_idx%rs(i))%sf(j, k, l) = muR
                    q_prim_vf(bub_idx%vs(i))%sf(j, k, l) = muV
                    if (.not. polytropic) then
                        q_prim_vf(bub_idx%ps(i))%sf(j, k, l) = patch_icpp(patch_id)%p0
                        q_prim_vf(bub_idx%ms(i))%sf(j, k, l) = patch_icpp(patch_id)%m0
                    end if
                end if
            end do
        end if

        ! Density and the specific heat ratio and liquid stiffness functions
        ! call s_convert_species_to_mixture_variables( &
        ! print*, 'third call to convert to mixture'
        call s_convert_to_mixture_variables( &
            q_prim_vf, j, k, l, &
            patch_icpp(smooth_patch_id)%rho, &
            patch_icpp(smooth_patch_id)%gamma, &
            patch_icpp(smooth_patch_id)%pi_inf)

        ! ==================================================================

        ! Pressure
        q_prim_vf(E_idx)%sf(j, k, l) = &
            (eta*patch_icpp(patch_id)%pres &
             + (1d0 - eta)*orig_prim_vf(E_idx))

        ! Volume fractions \alpha
        do i = adv_idx%beg, adv_idx%end
            q_prim_vf(i)%sf(j, k, l) = &
                eta*patch_icpp(patch_id)%alpha(i - E_idx) &
                + (1d0 - eta)*orig_prim_vf(i)
        end do

        ! Elastic Shear Stress
        if (hypoelasticity) then
            do i = 1, (stress_idx%end - stress_idx%beg) + 1
                q_prim_vf(i + stress_idx%beg - 1)%sf(j, k, l) = &
                    (eta*patch_icpp(patch_id)%tau_e(i) &
                     + (1d0 - eta)*orig_prim_vf(i + stress_idx%beg - 1))
            end do
        end if

        if (mpp_lim .and. bubbles) then
            !adjust volume fractions, according to modeled gas void fraction
            alf_sum%sf = 0d0
            do i = adv_idx%beg, adv_idx%end - 1
                alf_sum%sf = alf_sum%sf + q_prim_vf(i)%sf
            end do

            do i = adv_idx%beg, adv_idx%end - 1
                q_prim_vf(i)%sf = q_prim_vf(i)%sf*(1.d0 - q_prim_vf(alf_idx)%sf) &
                                  /alf_sum%sf
            end do
        end if

        ! Partial densities \alpha \rho
        if (model_eqns /= 4) then
            !mixture density is an input
            do i = 1, cont_idx%end
                q_prim_vf(i)%sf(j, k, l) = &
                    eta*patch_icpp(patch_id)%alpha_rho(i) &
                    + (1d0 - eta)*orig_prim_vf(i)
            end do
        else
            !get mixture density from pressure via Tait EOS
            pi_inf = fluid_pp(1)%pi_inf
            gamma = fluid_pp(1)%gamma
            lit_gamma = (1.d0 + gamma)/gamma

            ! \rho = (( p_l + pi_inf)/( p_ref + pi_inf))**(1/little_gam) * rhoref(1-alf)
            q_prim_vf(1)%sf(j, k, l) = &
                (((q_prim_vf(E_idx)%sf(j, k, l) + pi_inf)/(pref + pi_inf))**(1/lit_gamma))* &
                rhoref*(1 - q_prim_vf(alf_idx)%sf(j, k, l))
        end if

        ! Density and the specific heat ratio and liquid stiffness functions
        ! call s_convert_species_to_mixture_variables(q_prim_vf, j, k, l, &
        ! print*, 'fourth call to convert to mixture'
        call s_convert_to_mixture_variables(q_prim_vf, j, k, l, &
                                                    rho, gamma, pi_inf)

        ! Velocity
        do i = 1, E_idx - mom_idx%beg
            q_prim_vf(i + cont_idx%end)%sf(j, k, l) = &
                (eta*patch_icpp(patch_id)%vel(i) &
                 + (1d0 - eta)*orig_prim_vf(i + cont_idx%end))
        end do

        ! Set streamwise velocity to hypertangent function of y
         if (vel_profile) then
             q_prim_vf(1 + cont_idx%end)%sf(j, k, l) = &
                 (eta*patch_icpp(patch_id)%vel(1)*tanh(y_cc(k)) &
                 + (1d0 - eta)*orig_prim_vf(1 + cont_idx%end))
         end if
   
        ! Set partial pressures to mixture pressure for the 6-eqn model
        if (model_eqns == 3) then
            do i = internalEnergies_idx%beg, internalEnergies_idx%end
                q_prim_vf(i)%sf(j, k, l) = q_prim_vf(E_idx)%sf(j, k, l)
            end do
        end if

        ! Smoothed bubble variables
        if (bubbles) then
            do i = 1, nb
                muR = R0(i)*patch_icpp(patch_id)%r0 ! = 1*R0(i)
                muV = V0(i)*patch_icpp(patch_id)%v0 ! = 1*V0(i)
                if (qbmm) then
                    ! Initialize the moment set
                    if (dist_type == 1) then
                        q_prim_vf(bub_idx%fullmom(i, 0, 0))%sf(j, k, l) = 1d0
                        q_prim_vf(bub_idx%fullmom(i, 1, 0))%sf(j, k, l) = muR
                        q_prim_vf(bub_idx%fullmom(i, 0, 1))%sf(j, k, l) = muV
                        q_prim_vf(bub_idx%fullmom(i, 2, 0))%sf(j, k, l) = muR**2 + sigR**2
                        q_prim_vf(bub_idx%fullmom(i, 1, 1))%sf(j, k, l) = muR*muV + rhoRV*sigR*sigV
                        q_prim_vf(bub_idx%fullmom(i, 0, 2))%sf(j, k, l) = muV**2 + sigV**2
                    else if (dist_type == 2) then
                        q_prim_vf(bub_idx%fullmom(i, 0, 0))%sf(j, k, l) = 1d0
                        q_prim_vf(bub_idx%fullmom(i, 1, 0))%sf(j, k, l) = dexp((sigR**2)/2d0)*muR
                        q_prim_vf(bub_idx%fullmom(i, 0, 1))%sf(j, k, l) = muV
                        q_prim_vf(bub_idx%fullmom(i, 2, 0))%sf(j, k, l) = dexp((sigR**2)*2d0)*(muR**2)
                        q_prim_vf(bub_idx%fullmom(i, 1, 1))%sf(j, k, l) = dexp((sigR**2)/2d0)*muR*muV
                        q_prim_vf(bub_idx%fullmom(i, 0, 2))%sf(j, k, l) = muV**2 + sigV**2
                    end if
                else
                    ! q_prim_vf(bub_idx%rs(i))%sf(j,k,l) = &
                    !     (eta * R0(i)*patch_icpp(patch_id)%r0 &
                    !     + (1d0-eta)*orig_prim_vf(bub_idx%rs(i)))
                    ! q_prim_vf(bub_idx%vs(i))%sf(j,k,l) = &
                    !     (eta * V0(i)*patch_icpp(patch_id)%v0 &
                    !     + (1d0-eta)*orig_prim_vf(bub_idx%vs(i)))
                    q_prim_vf(bub_idx%rs(i))%sf(j, k, l) = muR
                    q_prim_vf(bub_idx%vs(i))%sf(j, k, l) = muV

                    if (.not. polytropic) then
                        q_prim_vf(bub_idx%ps(i))%sf(j, k, l) = patch_icpp(patch_id)%p0
                        q_prim_vf(bub_idx%ms(i))%sf(j, k, l) = patch_icpp(patch_id)%m0
                    end if

                end if
            end do
        end if

        if (mpp_lim .and. bubbles) then
            !adjust volume fractions, according to modeled gas void fraction
            alf_sum%sf = 0d0
            do i = adv_idx%beg, adv_idx%end - 1
                alf_sum%sf = alf_sum%sf + q_prim_vf(i)%sf
            end do

            do i = adv_idx%beg, adv_idx%end - 1
                q_prim_vf(i)%sf = q_prim_vf(i)%sf*(1.d0 - q_prim_vf(alf_idx)%sf) &
                                  /alf_sum%sf
            end do
        end if

        if (bubbles .and. (.not. polytropic) .and. (.not. qbmm)) then
            do i = 1, nb
                if (q_prim_vf(bub_idx%ps(i))%sf(j, k, l) == dflt_real) then
                    q_prim_vf(bub_idx%ps(i))%sf(j, k, l) = pb0(i)
                    print *, 'setting to pb0'
                end if
                if (q_prim_vf(bub_idx%ms(i))%sf(j, k, l) == dflt_real) then
                    q_prim_vf(bub_idx%ms(i))%sf(j, k, l) = mass_v0(i)
                end if
            end do
        end if
      
        ! Updating the patch identities bookkeeping variable
        if (1d0 - eta < 1d-16) patch_id_fp(j, k, l) = patch_id

    end subroutine s_assign_patch_species_primitive_variables! ------------

    subroutine s_finialize_assign_variables_module

        ! Nullifying procedure pointer to the subroutine assigning either
        ! the patch mixture or species primitive variables to a cell in the
        ! computational domain
        s_assign_patch_primitive_variables => null()

    end subroutine s_finialize_assign_variables_module

end module
