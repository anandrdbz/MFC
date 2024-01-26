!>
!! @file m_fftw.f90
!! @brief Contains module m_fftw

#:include 'macros.fpp'

!> @brief The module contains the subroutines for the FFT routines
module m_fftw

    ! Dependencies =============================================================
    use, intrinsic :: iso_c_binding

    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Message passing interface (MPI) module proxy

    ! ==========================================================================

    implicit none

    private; public :: s_initialize_fftw_module, &
 s_apply_fourier_filter, &
 s_finalize_fftw_module



contains

    !>  The purpose of this subroutine is to create the fftw plan
        !!      that will be used in the forward and backward DFTs when
        !!      applying the Fourier filter in the azimuthal direction.
    subroutine s_initialize_fftw_module() ! ----------------------------------


    end subroutine s_initialize_fftw_module ! ------------------------------

    !>  The purpose of this subroutine is to apply a Fourier low-
        !!      pass filter to the flow variables in the azimuthal direction
        !!      to remove the high-frequency content. This alleviates the
        !!      restrictive CFL condition arising from cells near the axis.
    subroutine s_apply_fourier_filter(q_cons_vf) ! --------------------------

        type(scalar_field), dimension(sys_size), intent(INOUT) :: q_cons_vf

        integer :: Nfq !< Number of kept modes

        integer :: i, j, k, l !< Generic loop iterators

        ! Restrict filter to processors that have cells adjacent to axis
    end subroutine s_apply_fourier_filter ! --------------------------------

    !>  The purpose of this subroutine is to destroy the fftw plan
        !!      that will be used in the forward and backward DFTs when
        !!      applying the Fourier filter in the azimuthal direction.
    subroutine s_finalize_fftw_module() ! ------------------------------------


    end subroutine s_finalize_fftw_module ! --------------------------------

end module
