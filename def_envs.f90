        module fit_snap_class
        use lists_class
        use proj_disp_class
        implicit none

        TYPE data_file
         integer                        :: frames
         integer                        :: nats
         character (len=100)            :: inp_data
         character (len=100)            :: inp_tens
         character (len=100)            :: inp_traj
         complex(8), allocatable        :: tens(:,:)
         character (len=5), allocatable :: label(:)
         double precision, allocatable  :: x(:,:)
         double precision, allocatable  :: x0(:)
         double precision, allocatable  :: mass(:)
         type(list), allocatable        :: neigh(:)
         type(molecule), allocatable    :: envs(:)
         double precision               :: weight
         integer, allocatable           :: cls(:)
         contains
         procedure                      :: set_neigh
         procedure                      :: set_envs
         procedure                      :: get_envs
        end TYPE data_file

        TYPE system
         type(data_file), allocatable   :: data(:)
         integer                        :: ndata
         character (len=100)            :: inp
         character (len=100)            :: inp_fix
         character (len=100)            :: inp_fit
         integer                        :: tot_frames
         integer                        :: npar2fit
         contains
         procedure                      :: read_sys
        END TYPE system

        TYPE kernel_kind
         integer                        :: nenvs=0
         double precision               :: sigma=0
         double precision, allocatable  :: B(:,:)
         double precision, allocatable  :: x0(:,:)
         contains
         procedure                      :: dealloc => dealloc_kernel_kind
        END TYPE kernel_kind

        TYPE kernel_global
         integer                        :: nkinds=0
         type(kernel_kind), allocatable :: K(:)
         contains
         procedure                      :: dealloc => dealloc_kernel_global
        END TYPE kernel_global
        
        contains

        subroutine dealloc_kernel_global(this)
        implicit none
        class(kernel_global)    :: this
        integer                 :: i
         if( allocated(this%K) ) then 
          do i=1,this%nkinds
           call this%K(i)%dealloc()
          enddo
          deallocate(this%K)
         endif
         this%nkinds=0
        return
        end subroutine dealloc_kernel_global

        subroutine dealloc_kernel_kind(this)
        implicit none
        class(kernel_kind)      :: this
         if( allocated(this%B) ) deallocate(this%B)
         this%nenvs=0
         this%sigma=0
        return
        end subroutine dealloc_kernel_kind

        subroutine read_sys(sys,data_file,tens_order)
        use rotations_class
        implicit none
        class(system)           :: sys
        integer                 :: i,l,v,k,m,n,tens_order,Tdim,Tdim2
        character(len=100)      :: label,data_file
        double precision, allocatable :: coeff(:,:)
        complex(8), allocatable   :: wig(:,:)
        double precision :: rot(3,3),vec(3),x(3),beta
        integer :: s1


         if(allocated(sys%data))then
          do i=1,sys%ndata 
           if(allocated(sys%data(i)%x)) deallocate(sys%data(i)%x)
           if(allocated(sys%data(i)%tens)) deallocate(sys%data(i)%tens)
          enddo 
          deallocate(sys%data)
         endif   

         allocate(coeff(2,2*tens_order+1))

         open(8,file=trim(data_file))

         sys%tot_frames=0
         sys%npar2fit=0

         read(8,*) sys%ndata    

         allocate(sys%data(sys%ndata))

         do i=1,sys%ndata

          read(8,*) sys%data(i)%inp_data,sys%data(i)%frames,sys%data(i)%inp_traj,sys%data(i)%inp_tens&
                    ,sys%data(i)%weight

          open(12,file=sys%data(i)%inp_traj)

          allocate(sys%data(i)%tens(sys%data(i)%frames,(2*tens_order+1)))
          open(13,file=sys%data(i)%inp_tens)

          do l=1,sys%data(i)%frames

           sys%tot_frames=sys%tot_frames+1

           read(12,*) sys%data(i)%nats
           read(12,*) 

           if(.not.allocated(sys%data(i)%label))then
            allocate(sys%data(i)%label(sys%data(i)%nats))
           endif
           if(.not.allocated(sys%data(i)%x))then
            allocate(sys%data(i)%x(sys%data(i)%frames,3*sys%data(i)%nats))
           endif

!           vec=0.0d0
!           vec(2)=1.0d0
!           call random_number(beta)
!           beta=beta*90.0d0-45.0d0
!           call Rot_Cart(beta*acos(-1.0d0)/180.0d0,vec,rot)
!           call Rot_Wig(dble(tens_order),0.0d0,beta*acos(-1.0d0)/180.0d0,0.0d0,wig)

           v=1
           do k=1,sys%data(i)%nats             
            read(12,*) sys%data(i)%label(k),sys%data(i)%x(l,v),sys%data(i)%x(l,v+1),sys%data(i)%x(l,v+2)             
!            x=0.0d0
!            do s1=1,3
!             x(s1)=x(s1)+rot(s1,1)*sys%data(i)%x(l,v)+&
!                         rot(s1,2)*sys%data(i)%x(l,v+1)+&
!                         rot(s1,3)*sys%data(i)%x(l,v+2)
!            enddo
!            sys%data(i)%x(l,v)=x(1)
!            sys%data(i)%x(l,v+1)=x(2)
!            sys%data(i)%x(l,v+2)=x(3)
            v=v+3
           enddo

           read(13,*) (coeff(1,Tdim),coeff(2,Tdim),Tdim=1,(2*tens_order+1))

!           do Tdim=1,2*tens_order+1
!            sys%data(i)%tens(l,Tdim)=(0.0d0,0.0d0)
!            do Tdim2=1,2*tens_order+1
!             sys%data(i)%tens(l,Tdim)=sys%data(i)%tens(l,Tdim)+conjg(wig(Tdim,Tdim2))*cmplx(coeff(1,Tdim2),coeff(2,Tdim2),8)
!            enddo
!           enddo

           do Tdim=1,2*tens_order+1
            sys%data(i)%tens(l,Tdim)=cmplx(coeff(1,Tdim),coeff(2,Tdim),8)
           enddo

           sys%npar2fit=sys%npar2fit+(2*tens_order+1)

          enddo

          close(12)
          close(13)

         enddo

         close(8)
         deallocate(coeff)

        return
        end subroutine read_sys

        subroutine set_neigh(sys,gen_cutoff)
        implicit none
        class(data_file)              :: sys
        integer                       :: s,l,i,v,v1,v2,k,t
        double precision              :: dist,gen_cutoff

         if (allocated(sys%neigh)) then
          do t=1,size(sys%neigh)
           call sys%neigh(t)%delete()
          enddo
          deallocate(sys%neigh)
         endif
         allocate(sys%neigh(sys%nats))

         do t=1,sys%nats        
          call sys%neigh(t)%init()
          do l=1,sys%nats
           dist=0.0d0
           do s=1,3
            v1=(t-1)*3+s
            v2=(l-1)*3+s
            dist=dist+(sys%x0(v1)-sys%x0(v2))**2
           enddo
           dist=sqrt(dist)
           if(dist.lt.gen_cutoff)then
            call sys%neigh(t)%add_node(l)
           endif
          enddo
          call sys%neigh(t)%reboot()
         enddo

        return
        end subroutine set_neigh

        subroutine set_envs(sys)
        implicit none
        class(data_file)              :: sys
        integer                       :: gen_cutoff,s,l,i,v,v1,v2,k,t
        double precision              :: dist
        double precision, allocatable :: env(:,:),mass(:)

         if (allocated(sys%envs)) deallocate(sys%envs)
         if (allocated(sys%mass)) deallocate(sys%mass)
         allocate(sys%envs(sys%nats))
         allocate(sys%mass(sys%nats))
         sys%mass=1.0d0

         do t=1,sys%nats        
          call sys%neigh(t)%reboot()
          if (allocated(env)) deallocate(env)
          if (allocated(mass)) deallocate(mass)
          allocate(env(sys%neigh(t)%nelem,3))
          allocate(mass(sys%neigh(t)%nelem))
          do l=1,sys%neigh(t)%nelem                   
           call sys%neigh(t)%rd_val(k)
           mass(l)=sys%mass(k)
           do s=1,3
            v2=(k-1)*3+s
            env(l,s)=sys%x0(v2)
           enddo
           call sys%neigh(t)%skip()
          enddo
          call sys%neigh(t)%reboot()
          call sys%envs(t)%def_mol(env,mass)
         enddo
  
         if (allocated(env)) deallocate(env)

        return
        end subroutine set_envs

        subroutine get_envs(sys,i)
        implicit none
        class(data_file)              :: sys
        integer                       :: gen_cutoff,s,l,i,v,v1,v2,k,t
        double precision              :: dist
        double precision, allocatable :: env(:,:)

!         do t=1,sys%nats        
          call sys%neigh(1)%reboot()
          if (allocated(env)) deallocate(env)
          allocate(env(sys%neigh(1)%nelem,3))
          do l=1,sys%neigh(1)%nelem                   
           call sys%neigh(1)%rd_val(k)
           do s=1,3
            v2=(k-1)*3+s
            env(l,s)=sys%x(i,v2)
           enddo
           call sys%neigh(1)%skip()
          enddo
          call sys%neigh(1)%reboot()
          call sys%envs(1)%def_mol_dist(env)
          call sys%envs(1)%proj_disp()
          call sys%envs(1)%get_euler()              
!          write(*,*) sys%envs(t)%alpha,sys%envs(t)%beta,sys%envs(t)%gamma
!          write(*,*) sys%neigh(t)%nelem
!          write(*,*)
!          do l=1,sys%neigh(t)%nelem
!           write(*,*) 'H',sys%envs(t)%cart_eq(l,1),sys%envs(t)%cart_eq(l,2),sys%envs(t)%cart_eq(l,3)
!          enddo
!          write(*,*) sys%neigh(t)%nelem
!          write(*,*)
!          do l=1,sys%neigh(t)%nelem
!           write(*,*) 'H',env(l,1),env(l,2),env(l,3)
!          enddo
!          write(*,*) sys%neigh(t)%nelem
!          write(*,*)
!          do l=1,sys%neigh(t)%nelem
!           write(*,*) 'H',sys%envs(t)%cart_rot(l,1),sys%envs(t)%cart_rot(l,2),sys%envs(t)%cart_rot(l,3)
!          enddo
!          stop
!         enddo
  
         if (allocated(env)) deallocate(env)

        return
        end subroutine get_envs

        end module fit_snap_class
