#:def LOG(expr)
#ifdef MFC_DEBUG
    block
        use iso_fortran_env, only: output_unit

        print *, '${_FILE_.split('/')[-1]}$:${_LINE_}$: ', ${expr}$
        call flush(output_unit)
    end block
#endif
#:enddef

#:def ALLOCATE(*args)
    allocate(${', '.join(args)}$)
    !$acc enter data create(${', '.join([ arg.split('(')[0] for arg in args ])}$)
#:enddef ALLOCATE

#:def DEALLOCATE(*args)
    deallocate(${', '.join(args)}$)
    !$acc exit data delete(${', '.join(args)}$)
#:enddef DEALLOCATE

#define t_vec3   real(kind(0d0)), dimension(1:3)
#define t_mat4x4 real(kind(0d0)), dimension(1:4,1:4)
