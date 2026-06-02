#:mute
#:def DTIMES(text)
${ ', '.join([text] * NDIM) }$
#:enddef

#:def DINDEX(array)
#:if NDIM == 2
${array}$(1), ${array}$(2)
#:elif NDIM == 3
${array}$(1), ${array}$(2), ${array}$(3)
#:endif
#:enddef

#:def KJI_LOOP_1_to_array(array)
#:if NDIM == 2
j = 1, ${array}$(2); do i = 1, ${array}$(1)
#:elif NDIM == 3
k = 1, ${array}$(3); do j = 1, ${array}$(2); do i = 1, ${array}$(1)
#:endif
#:enddef

#:def KJI_LOOP_array_to_array(lo, hi)
#:if NDIM == 2
j = ${lo}$(2), ${hi}$(2); do i = ${lo}$(1), ${hi}$(1)
#:elif NDIM == 3
k = ${lo}$(3), ${hi}$(3); do j = ${lo}$(2), ${hi}$(2); do i = ${lo}$(1), ${hi}$(1)
#:endif
#:enddef

#:if NDIM == 2
#:set KJI_CLOSE_LOOP='end do'
#:elif NDIM == 3
#:set KJI_CLOSE_LOOP='end do; end do'
#:endif

#:if NDIM == 2
#:set IJK='i, j'
#:set XYZ='x, y'
#:elif NDIM == 3
#:set IJK='i, j, k'
#:set XYZ='x, y, z'
#:endif

#:endmute
