        subroutine get_clusters(nclusters)
        use common_var
        use LAMMPS
        use, intrinsic :: ISO_C_binding, only : C_double, C_ptr, C_int
        implicit none
        double precision                 :: chi_val_ener,small,dist1,summ
        complex(8)                       :: ener,dev
        integer                          :: l,i,k,n,v,nclusters,j,m
        integer                          :: npar,nkinds,iter
        real (C_double), pointer         :: bisp(:,:) => null()
        real (c_double), pointer         :: kind_nat(:) => null()
        real (C_double), pointer         :: id_dbl(:)=> null()
        integer, allocatable             :: map(:),id(:)
        character(len=1000)              :: word
        double precision, allocatable    :: B(:,:),Mean(:,:),dist(:),new_Mean(:,:)
        integer, allocatable             :: counter(:)
        logical                          :: check

         npar=55
         allocate(B(npar,sys%tot_frames))       
         allocate(Mean(npar,nclusters))
         allocate(new_Mean(npar,nclusters))
         allocate(dist(nclusters))
         allocate(counter(nclusters))

         n=1

         do i=1,sys%ndata

          allocate(sys%data(i)%cls(sys%data(i)%frames))

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

           summ=0.0d0
           do v=1,npar
            B(v,n)=bisp(v,map(1))
            summ=summ+B(v,n)**2
           enddo

           B(:,n)=B(:,n)/sqrt(summ)

           n=n+1

          enddo
         enddo

       ! Do 2D map of Similarity

         open(121,file='kernel.dat')

         n=1
         do i=1,sys%ndata
          do l=1,sys%data(i)%frames
           k=1
           do j=1,sys%ndata
            do m=1,sys%data(j)%frames
             
             if(k.le.n)then

             dist1=0.0d0 
             do v=1,npar
              dist1=dist1+(B(v,n)-B(v,k))**2
             enddo

             write(121,*) n,k,dist1
        
             endif

            k=k+1
            enddo
           enddo
           write(121,*) 
           n=n+1
          enddo
         enddo
         close(121)

       ! Assign starting means

         Mean(:,1)=B(:,1)
         Mean(:,2)=B(:,2)
         check=.false.
         iter=1

         do while (iter.lt.40 .or. check)

          n=1
          do i=1,sys%ndata
           do l=1,sys%data(i)%frames
 
            do k=1,nclusters
             dist(k)=0.0d0
             do v=1,npar
              dist(k)=dist(k)+(B(v,n)-Mean(v,k))**2
             enddo
            enddo

        ! Assign elements to clusters

            small=dist(1)
            sys%data(i)%cls(l)=1
            do k=2,nclusters
             if(dist(k).lt.small)then
              small=dist(k)
              sys%data(i)%cls(l)=k
             endif
            enddo
            write(*,*) i,l,sys%data(i)%cls(l)

            n=n+1
           enddo
          enddo

        ! update means and cycle

          new_Mean=0.0d0
          counter=0

          n=1

          do i=1,sys%ndata
           do l=1,sys%data(i)%frames
 
            k=sys%data(i)%cls(l)
            counter(k)=counter(k)+1
            do v=1,npar
             new_Mean(v,k)=new_Mean(v,k)+B(v,n)
            enddo

            n=n+1
           enddo
          enddo

          do k=1,nclusters
           new_Mean(:,k)=new_Mean(:,k)/counter(k)
          enddo         
           
          Mean=new_Mean
          iter=iter+1
         enddo

        return
        end subroutine get_clusters
