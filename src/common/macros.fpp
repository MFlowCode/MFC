#:def LOG(expr)
#ifdef MFC_DEBUG
print *, '${_FILE_.split('/')[-1]}$:${_LINE_}$', ${expr}$
#endif
#:enddef


#:def ALLOCATE(*args)

! ================================================
! ==== BEGIN ALLOCATE (src/common/macros.fpp) ====
! ================================================

#! @:LOG({'@:ALLOCATE(${', '.join(args)}$)'})

allocate(${', '.join(args)}$)

!$acc enter data create(${', '.join(args)}$)

! ================================================
! ====  END  ALLOCATE (src/common/macros.fpp) ====
! ================================================

#:enddef ALLOCATE



#:def DEALLOCATE(*args)

! ==================================================
! ==== BEGIN DEALLOCATE (src/common/macros.fpp) ====
! ==================================================

#! @:LOG({'@:DEALLOCATE(${', '.join(args)}$)'})

deallocate(${', '.join(args)}$)

!$acc exit data delete(${', '.join(args)}$)

! ==================================================
! ====  END  DEALLOCATE (src/common/macros.fpp) ====
! ==================================================

#:enddef DEALLOCATE
