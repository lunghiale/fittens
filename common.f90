        module common_var
        use, intrinsic :: ISO_C_binding, only : C_double, C_ptr, C_int
        use fit_snap_class

         TYPE(system)                   :: sys
         type (C_ptr), allocatable      :: lmp(:)
         type(kernel_global)            :: kernel
         logical                        :: cs=.false.,skip_fit=.false.,sph=.true.
         double precision               :: cm_val=1.0d0,thr_kernel=0.5d0
         integer                        :: tens_order=0
                                                
        end module common_var

