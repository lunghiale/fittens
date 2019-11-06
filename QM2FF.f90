        program fitsnap
        use fit_snap_class
        use common_var
        implicit none
        integer                        :: l
        character (len=100)            :: command,input,datas,output,new_datas
       
         if(iargc().eq.0)then
          write(*,*) 'FitTens Usage:'
          write(*,*) '-datas    : Reference data'
          write(*,*) '-inp      : Lammps input file'
          write(*,*) '-pot_fit  : SNAP potential to fit' 
          write(*,*) '-pot_fix  : Static potential'
          write(*,*) '-tensor   : order of the tensor'
          write(*,*) '-compress : L2 regularization parameter'
         stop
         endif   

         do l=1,iargc()

          call getarg(l,command)

          if(trim(command).eq.'-datas')then        
           call getarg(l+1,command)
           read(command,*) datas
          endif
          if(trim(command).eq.'-inp')then        
           call getarg(l+1,command)
           read(command,*) sys%inp
          endif
          if(trim(command).eq.'-pot_fit')then        
           call getarg(l+1,command)
           read(command,*) sys%inp_fit
          endif
          if(trim(command).eq.'-pot_fix')then        
           call getarg(l+1,command)
           read(command,*) sys%inp_fix
          endif
          if(trim(command).eq.'-tensor')then        
           call getarg(l+1,command)
           read(command,*) tens_order
          endif
          if(trim(command).eq.'-skip_fit')then        
           skip_fit=.true.
          endif
          if(trim(command).eq.'-compress')then        
           cs=.true.
           call getarg(l+1,command)
           read(command,*) cm_val
          endif

         enddo

         call sys%read_sys(datas,tens_order)
         call init_lammps_lib
!         call get_clusters(2)
         call get_lsmf_snap
         call get_chi2

        return
        end program fitsnap


        subroutine get_chi2 
        use fit_snap_class
        use rotations_class
        use proj_disp_class
        use common_var
        use LAMMPS
        use, intrinsic :: ISO_C_binding, only : C_double, C_ptr, C_int
        implicit none
        double precision                 :: chi_val_ener,avg_sph(5)
        complex(8)                       :: ener,dev,fi
        integer                          :: l,i,k,n,v,Tdim,Tdim2,s
        integer                          :: npar,nkinds
        real (C_double), pointer         :: bisp(:,:) => null()
        real (c_double), pointer         :: kind_nat(:) => null()
        real (C_double), pointer         :: id_dbl(:)=> null()  
        real (C_double), pointer         :: f(:,:) => null()
        integer, allocatable             :: map(:),id(:)
        character(len=1000)              :: word
        double precision, allocatable    :: bcoeff(:,:)
        complex(8), allocatable          :: wig(:,:),B(:,:,:)
       
         open(411,file='snapcoeff_ref')
         read(411,*) nkinds,npar,tens_order

         allocate(B(npar,nkinds,2*tens_order+1))
         allocate(bcoeff(2,2*tens_order+1))
         do i=1,nkinds 
          read(411,*)
          do l=1,npar
           read(411,*) (bcoeff(1,Tdim),bcoeff(2,Tdim),Tdim=1,(2*tens_order+1))
           do Tdim=1,2*tens_order+1
            B(l,i,Tdim)=cmplx(bcoeff(1,Tdim),bcoeff(2,Tdim),8)
           enddo
          enddo
         enddo
         deallocate(bcoeff)
         close(411)

!!      calcola chi2
         do Tdim=1,(2*tens_order+1)

          chi_val_ener=0.0d0
          if(Tdim.lt.10) write(word,'(I1)') Tdim
          if(Tdim.lt.100 .and. Tdim.ge.10) write(word,'(I2)') Tdim
          open(222,file='tens_'//trim(word)//'_rmse.dat')
          if(sph) open(333,file='tens_'//trim(word)//'_sph.dat')

          do i=1,sys%ndata
           do l=1,sys%data(i)%frames

            call lammps_scatter_atoms (lmp(i),'x',sys%data(i)%x(l,1:3*sys%data(i)%nats))
            call lammps_command (lmp(i), 'run 0')
            !!!!!
            call lammps_extract_compute (kind_nat, lmp(i), 'type', 1, 1)
            call lammps_extract_compute (bisp, lmp(i), 'sna_e', 1, 2)
            if(.not. allocated(id)) allocate(id(sys%data(i)%nats))
            if(.not. allocated(map)) allocate(map(sys%data(i)%nats))
            call lammps_extract_compute (id_dbl, lmp(i), 'id', 1, 1)
          
            id=INT(id_dbl)
            id_dbl=>null()

            do k=1,sys%data(i)%nats          
             map(id(k))=k
            enddo

            ener=(0.0d0,0.0d0)

            call sys%data(i)%get_envs(l)
            call Rot_Wig(dble(tens_order),sys%data(i)%envs(1)%alpha,sys%data(i)%envs(1)%beta,sys%data(i)%envs(1)%gamma,wig)
!            call Rot_Wig(dble(tens_order),0.0d0,0.0d0,0.0d0,wig)

            do k=1,sys%data(i)%nats         
             do Tdim2=1,(2*tens_order+1)
              ener=ener+conjg(wig(Tdim,Tdim2))*B(1,nint(kind_nat(k)),Tdim2)
             enddo
             do v=2,size(B,1)
              do Tdim2=1,(2*tens_order+1)
               ener=ener+conjg(wig(Tdim,Tdim2))*B(v,nint(kind_nat(k)),Tdim2)*bisp(v-1,map(k))
              enddo
             enddo

            enddo

            dev=ener-sys%data(i)%tens(l,Tdim)        
            chi_val_ener=chi_val_ener+dble(dev*conjg(dev))

            write(222,*) i,l,dble(ener),aimag(ener),dble(sys%data(i)%tens(l,Tdim)),&
                             aimag(sys%data(i)%tens(l,Tdim)),dble(dev),aimag(dev)


            if(sph)then

             avg_sph=0.0d0
             
             call lammps_extract_compute (f, lmp(i), 'sna_f', 1, 2)
             
             do k=1,sys%data(i)%nats         
              do n=1,nkinds
               s=(n-1)*(3*(npar-1))+1
               fi=(0.0d0,0.0d0)
               do v=2,size(B,1)               
                fi=fi+f(s,map(k))*B(v,nint(kind_nat(k)),Tdim)
                s=s+1
               enddo
               avg_sph(Tdim)=avg_sph(Tdim)+dble(fi*conjg(fi))
              enddo
             enddo

             do k=1,sys%data(i)%nats         
              do n=1,nkinds
               s=(n-1)*(3*(npar-1))+1+(npar-1)
               fi=(0.0d0,0.0d0)
               do v=2,size(B,1)               
                fi=fi+f(s,map(k))*B(v,nint(kind_nat(k)),Tdim)
                s=s+1
               enddo
               avg_sph(Tdim)=avg_sph(Tdim)+dble(fi*conjg(fi))
              enddo
             enddo

             do k=1,sys%data(i)%nats         
              do n=1,nkinds
               s=(n-1)*(3*(npar-1))+1+2*(npar-1)
               fi=(0.0d0,0.0d0)
               do v=2,size(B,1)               
                fi=fi+f(s,map(k))*B(v,nint(kind_nat(k)),Tdim)
                s=s+1
               enddo
               avg_sph(Tdim)=avg_sph(Tdim)+dble(fi*conjg(fi))
              enddo
             enddo

             f=>null()

             write(333,*) i,l,avg_sph(Tdim)/(sys%data(i)%nats*3)

            endif

           enddo      ! ciclo su frames
          enddo   ! ciclo su data

          write(222,*) 'RMSE:',sqrt(chi_val_ener/sys%tot_frames)
          close(222)
          write(*,*) Tdim,' RMSE:',sqrt(chi_val_ener/sys%tot_frames)
          if(sph) close(333)

         enddo

         do i=1,sys%ndata
          call lammps_close (lmp(i))
         enddo

         if(allocated(lmp)) deallocate(lmp)
         if(allocated(map)) deallocate(map)
         if(allocated(id)) deallocate(id)

        return
        end subroutine get_chi2


!        subroutine get_D(x0,D)
!        use fit_snap_class
!        use rotations_class
!        use proj_disp_class
!        use common_var
!        use LAMMPS
!        use, intrinsic :: ISO_C_binding, only : C_double, C_ptr, C_int
!        implicit none
!        double precision                 :: chi_val_ener,avg_sph(5)
!        complex(8)                       :: ener,dev,fi
!        integer                          :: l,i,k,n,v,Tdim,Tdim2,s
!        integer                          :: npar,nkinds
!        real (C_double), pointer         :: bisp(:,:) => null()
!        real (c_double), pointer         :: kind_nat(:) => null()
!        real (C_double), pointer         :: id_dbl(:)=> null()  
!        real (C_double), pointer         :: f(:,:) => null()
!        integer, allocatable             :: map(:),id(:)
!        character(len=1000)              :: word
!        double precision, allocatable    :: bcoeff(:,:)
!        complex(8), allocatable          :: wig(:,:),B(:,:,:),T2(:)
       
!         open(411,file='snapcoeff_ref')
!         read(411,*) nkinds,npar,tens_order

!         allocate(T2(2*tens_order+1))
!         allocate(B(npar,nkinds,2*tens_order+1))
!         allocate(bcoeff(2,2*tens_order+1))
!         do i=1,nkinds 
!          read(411,*)
!          do l=1,npar
!           read(411,*) (bcoeff(1,Tdim),bcoeff(2,Tdim),Tdim=1,(2*tens_order+1))
!           do Tdim=1,2*tens_order+1
!            B(l,i,Tdim)=cmplx(bcoeff(1,Tdim),bcoeff(2,Tdim),8)
!           enddo
!          enddo
!         enddo
!         deallocate(bcoeff)
!         close(411)

!         call lammps_scatter_atoms (lmp(1),'x',x0(l,1:3*sys%data(1)%nats))
!         call lammps_command (lmp(1), 'run 0')
!         call lammps_extract_compute (kind_nat, lmp(1), 'type', 1, 1)
!         call lammps_extract_compute (bisp, lmp(1), 'sna_e', 1, 2)
!         if(.not. allocated(id)) allocate(id(sys%data(1)%nats))
!         if(.not. allocated(map)) allocate(map(sys%data(1)%nats))
!         call lammps_extract_compute (id_dbl, lmp(1), 'id', 1, 1)
          
!         id=INT(id_dbl)
!         id_dbl=>null()

!         do k=1,sys%data(1)%nats          
!          map(id(k))=k
!         enddo

!         call sys%data(1)%get_envs(l)
!         call Rot_Wig(dble(tens_order),sys%data(i)%envs(1)%alpha,sys%data(1)%envs(1)%beta,sys%data(1)%envs(1)%gamma,wig)

!         do Tdim=1,(2*tens_order+1)

!          ener=(0.0d0,0.0d0)

!          do k=1,sys%data(1)%nats         
!           do Tdim2=1,(2*tens_order+1)
!            ener=ener+conjg(wig(Tdim,Tdim2))*B(1,nint(kind_nat(k)),Tdim2)
!           enddo
!           do v=2,size(B,1)
!            do Tdim2=1,(2*tens_order+1)
!             ener=ener+conjg(wig(Tdim,Tdim2))*B(v,nint(kind_nat(k)),Tdim2)*bisp(v-1,map(k))
!            enddo
!           enddo
!          enddo

!          T2(Tdim)=ener

!         enddo   ! ciclo Tdim         

!         call lammps_close (lmp(1))
!         if(allocated(lmp)) deallocate(lmp)
!         if(allocated(map)) deallocate(map)
!         if(allocated(id)) deallocate(id)

!        return
!        end subroutine get_D


