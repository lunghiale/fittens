        subroutine init_lammps_lib
        use fit_snap_class 
        use common_var
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
          write(*,*) label(i),type2fit(i),radii(i),cutoff(i),sigma(i)
         enddo         
         close(16)
         write(snap_string,*) gen_cutoff,'1.0000',bi_order,(radii(i),i=1,nkinds),(cutoff(i),i=1,nkinds),&
                 'quadraticflag ',quadflag
         write(snap_string2,*) (type2fit(i),i=1,nkinds)
         write(*,*) trim(snap_string)
         write(*,*) 'Fitting types: ',trim(snap_string2)

         open(333,file='snapparam')        
         write(333,*) 'rcutfac ',gen_cutoff
         write(333,*) 'twojmax ',bi_order
         write(333,*) 'quadraticflag ',quadflag
         write(333,*) 'rfac0 1.00000'
         write(333,*) 'rmin0 0'
         write(333,*) 'diagonalstyle 3'
         write(333,*) 'switchflag 1'

         close(333)

         open(222,file='snapcoeff')

         write(222,*) nkinds,npar
         do i =1,nkinds
          write(222,*) label(i),radii(i),cutoff(i)
          do n=1,npar
           write(222,*) 0.0000000
          enddo
         enddo
         flush(222)
         close(222)
        
         allocate(lmp(sys%ndata))

         do i=1,sys%ndata
           
          call lammps_open_no_mpi ('lmp -log none', lmp(i))
          call lammps_file (lmp(i),sys%inp)
          call lammps_command (lmp(i),'read_data '//trim(sys%data(i)%inp_data))
          call lammps_command (lmp(i),'group fitsnap type '//trim(snap_string2))
          call lammps_file (lmp(i),sys%inp_fix)
          call lammps_command (lmp(i), &
                 'compute sna_e all sna/atom '//trim(snap_string)//&
                 ' diagonal 3 rmin0 0 switchflag 1')
          if(sph) call lammps_command (lmp(i), &
                  'compute sna_f all snad/atom '//trim(snap_string)//&
                  ' diagonal 3 rmin0 0 switchflag 1')
          call lammps_command (lmp(i),&
               'compute type all property/atom type')
          call lammps_command (lmp(i),&
               'compute id all property/atom id')

         enddo

! Do kernel Neighbours list and  Environments

         call kernel%dealloc()

         kernel%nkinds=nkinds
         allocate(kernel%K(nkinds))
         do i=1,nkinds
          kernel%K(i)%sigma=sigma(i)
          kernel%K(i)%nenvs=0
         enddo

         do i=1,sys%ndata
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
          enddo

          call lammps_gather_atoms (lmp(i),'x',3,sys%data(i)%x0)
          call sys%data(i)%set_neigh(gen_cutoff)
          call sys%data(i)%set_envs()

          kind_count=0
          do t=1,sys%data(i)%nats           
           kind_count(nint(kind_nat(t)))=kind_count(nint(kind_nat(t)))+1
          enddo
          do t=1,nkinds
           kernel%K(t)%nenvs=kernel%K(t)%nenvs+kind_count(t)*sys%data(i)%frames
          enddo

         enddo

         do t=1,nkinds
          allocate(kernel%K(t)%B(kernel%K(t)%nenvs,size(bisp,1)))
         enddo

         do i=1,nkinds
          kernel%K(i)%nenvs=0
         enddo

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
           enddo

           do t=1,sys%data(i)%nats
            v=kernel%K(nint(kind_nat(map(t))))%nenvs+1
            do k=1,size(bisp,1)
             kernel%K(nint(kind_nat(map(t))))%B(v,k)=bisp(k,map(t))
            enddo
            kernel%K(nint(kind_nat(map(t))))%nenvs=&
                kernel%K(nint(kind_nat(map(t))))%nenvs+1
           enddo

           id_dbl=>null()
           bisp=>null()
           kind_nat=>null()

          enddo

          deallocate(id)
          deallocate(map)

         enddo
             
        return
        end subroutine init_lammps_lib
