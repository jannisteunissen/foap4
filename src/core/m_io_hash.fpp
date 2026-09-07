#:include 'definitions_ndim.fpp'
module m_io_hash_${NDIM}$d
  use iso_fortran_env
  implicit none

  type, public :: key_t
     integer :: x(${NDIM}$+1)
  end type key_t

#define FFH_KEY_TYPE type(key_t)
#define FFH_VAL_TYPE integer
#define FFH_CUSTOM_KEYS_EQUAL
#include "ffhash_inc.f90"

  pure logical function keys_equal(a, b)
    type(key_t), intent(in) :: a, b
    keys_equal = all(a%x == b%x)
  end function keys_equal

end module m_io_hash_${NDIM}$d
