        subroutine get_lsmf_snap
        use fit_snap_class 
        use lapack_diag_simm
        use lapack_inverse
        use common_var
        use rotations_class
        use lists_class
        use proj_disp_class
        use LAMMPS
        use, intrinsic :: ISO_C_binding, only : C_double, C_ptr, C_int
        implicit none
        integer                         :: l,i,k,n,v,s,t,v1,v2
        integer                         :: ss,inf,LWORK
        complex(8), allocatable         :: B(:),A(:,:),work(:)
        double precision                :: alpha,beta,gamma,error
        complex(8), allocatable         :: wig(:,:)

        real (C_double), pointer        :: bisp(:,:) => null()
        real (c_double), pointer        :: kind_nat(:) => null()
        real (C_double), pointer        :: id_dbl(:)=> null()
        integer, allocatable            :: map(:),id(:)

        integer                         :: nkinds,bi_order,npar,npar2fit,tot_frames
        integer                         :: nats2fit,quadflag,dimA,dimB,Tdim,Tdim2
        double precision                :: gen_cutoff
        character (len=10000)           :: snap_string,snap_string2,word
        double precision, allocatable   :: cutoff(:),radii(:),sigma(:)
        integer, allocatable            :: kind_count(:),type2fit(:)
        character (len=2),allocatable   :: label(:)

         tot_frames=sys%tot_frames         
         npar2fit=sys%npar2fit
      
         open(16,file=sys%inp_fit)
         read(16,*) gen_cutoff,bi_order,npar,quadflag
         read(16,*) nkinds
         allocate(label(nkinds))
         allocate(type2fit(nkinds))
         allocate(cutoff(nkinds))
         allocate(sigma(nkinds))
         allocate(radii(nkinds))
         allocate(kind_count(nkinds))
         kind_count=0
         do i=1,nkinds
          read(16,*) label(i),type2fit(i),radii(i),cutoff(i),sigma(i)
         enddo         
         close(16)

! Do fitting

         if(.not.skip_fit)then

          v=1
          do i=1,sys%ndata
           do l=1,sys%data(i)%frames

            call lammps_scatter_atoms (lmp(i),'x',sys%data(i)%x(l,1:3*sys%data(i)%nats))
            call lammps_command (lmp(i), 'run 0')
            call lammps_extract_compute (kind_nat, lmp(i), 'type', 1, 1)
            call lammps_extract_compute (bisp, lmp(i), 'sna_e', 1, 2)
            if(.not. allocated(id)) allocate(id(sys%data(i)%nats))
            if(.not. allocated(map)) allocate(map(sys%data(i)%nats))
            call lammps_extract_compute (id_dbl, lmp(i), 'id', 1, 1)
          
            id=INT(id_dbl)
            id_dbl=>null()

            do k=1,sys%data(i)%nats          
             map(id(k))=k
             writE(*,*) k,map(k),id(k)
            enddo

            if(l.eq.1)then

             npar=size(bisp,1)+1
             kind_count=0
            
             do t=1,sys%data(i)%nats
              kind_count(nint(kind_nat(t)))=kind_count(nint(kind_nat(t)))+1
             enddo

             if(.not.cs)then
              dimA=(2*tens_order+1)*npar*nkinds
              dimB=(npar2fit)*(2*tens_order+1)
             else
              dimA=(2*tens_order+1)*npar*nkinds
              dimB=(npar2fit+nkinds+(npar-1)*nkinds)*(2*tens_order+1)
             endif
             if(.not.allocated(B)) allocate(B(dimB))
              if(.not.allocated(A)) then
              allocate(A(dimB,dimA))
              A=(0.0d0,0.0d0)
             endif

            endif

            call sys%data(i)%get_envs(l)
            call Rot_Wig(dble(tens_order),sys%data(i)%envs(1)%alpha,sys%data(i)%envs(1)%beta,sys%data(i)%envs(1)%gamma,wig)
!            call Rot_Wig(dble(tens_order),0.0d0,0.0d0,0.0d0,wig)

            do Tdim=1,(2*tens_order+1)
             do t=1,sys%data(i)%nats

              do Tdim2=1,(2*tens_order+1)
               s=1+(nint(kind_nat(map(t)))-1)*npar+(Tdim2-1)*npar*nkinds
               A(v,s)=A(v,s)+conjg(wig(Tdim,Tdim2))
              enddo
                         
              do k=2,npar
               do Tdim2=1,(2*tens_order+1)
                s=k+(nint(kind_nat(t))-1)*npar+(Tdim2-1)*npar*nkinds
                A(v,s)=A(v,s)+bisp(k-1,map(t))*sys%data(i)%weight*conjg(wig(Tdim,Tdim2))
               enddo
              enddo

             enddo
             
             B(v)=sys%data(i)%tens(l,Tdim)*sys%data(i)%weight
             v=v+1
            
            enddo ! Tdim cycle        
            bisp=>null()
           enddo ! ciclo on frames

           if(allocated(map)) deallocate(map)
           if(allocated(id)) deallocate(id)
    
          enddo   ! ciclo su data


! compressive sensing

          if(cs)then

           do Tdim=1,(2*tens_order+1)
            s=1+(Tdim-1)*npar*nkinds
            do k=1,nkinds!-1
             B(v)=(0.0d0,0.0d0)
             A(v,s)=cmplx(cm_val,cm_val,8)
             s=s+npar
             v=v+1
            enddo
           enddo

           do Tdim=1,(2*tens_order+1)
            do k=1,nkinds
             do l=2,npar
              B(v)=(0.0d0,0.0d0)
              s=(k-1)*npar+l+(Tdim-1)*npar*nkinds
              A(v,s)=cmplx(cm_val,cm_val,8)
              v=v+1
             enddo
            enddo
           enddo

          endif

! Fitting 

          write(*,*) '######## Fitting'
          flush(6)
          lwork=dimB+64*dimB+1000
          allocate(work(lwork))
          call zgels('N',dimB,dimA,1,A,dimB,B,dimB,WORK,LWORK,inf)
          deallocate(work)
          if(inf.ne.0)then
           write(*,*) 'zgels failed',inf
           stop
          endif 
          write(*,*) '######## Done'
          error=0.0d0
          do i=dimA+1,dimB
           error=error+conjg(B(i))*B(i)
          enddo
          open(1111,file='rmse')
          write(1111,*) sqrt(error/dimB)
          flush(1111)
          close(1111)
          flush(6)

! Printing out         

          open(333,file='snapparam')        
          write(333,*) 'rcutfac ',gen_cutoff
          write(333,*) 'twojmax ',bi_order
          write(333,*) 'quadraticflag ',quadflag
          write(333,*) 'rfac0 1.00000'
          write(333,*) 'rmin0 0'
          write(333,*) 'diagonalstyle 3'
          write(333,*) 'switchflag 1'

          close(333)

          l=1
          open(222,file='snapcoeff_ref')
          write(222,*) nkinds,npar,tens_order
          do i =1,nkinds
           write(222,*) label(i),radii(i),cutoff(i)
           do n=1,npar
            write(222,*) (dble(B(l+(Tdim-1)*npar*nkinds)),&
                       aimag(B(l+(Tdim-1)*npar*nkinds)),Tdim=1,(2*tens_order+1))
             l=l+1
           enddo
          enddo
          close(222)

         endif ! skip_fit 
             
        return
        end subroutine get_lsmf_snap
