program hk2hr
!ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
!                                                                           s
!      The program converts H(k) (Hamiltonian in Reciprocal Space)          s
!              to H(r) (in  Direct Space) & computes hoppings               s
!                                                                           s
!           by Andrew Korneev, supervized by Alexander Poteryaev            s
!                                                                           s
!ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
!                                                                           s
!      Input:                                                               s
!      Hamiltonian file and system.am file                                  s
!                                                                           s
!      Output:                                                              s
!      hoppings.out                                                         s
!                                                                           s
!                     1 ---                                                 s
!          t_ij(r) = ___\     ik(R_i - R_j - transl_vec)                    s
!                     N /   e                            * H_ij(k)          s
!                       ---                                                 s
!                        k                                                  s
!ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss

  use hamiltonian_module
  use parameter_module
   

!============================================================================
!
!                       Variables initialization
!
!============================================================================

  implicit none

  integer :: iuham = 2, ierr                                  
  character(128) :: ham_name,system_name = 'system.am', scanner
  integer :: i, j, k, natoms, i_central, nblocks, dim, stopmark=0, stopmark1=0, nnbrs=0, ntransl, iatom, is, m1, m2, ikp
  real(8) :: brav_lattice(3,3), transl_vec(3), new_pos(3), b(3,3), M(3)
  real(8) :: alat, radius, distance, k_r, V, pi
  logical :: tag_exist
  real(8), allocatable :: at_position(:,:), dec_at_position(:,:), dec_k(:,:)
  complex(8), allocatable :: H_r(:,:,:)
  character(4), allocatable :: at_name(:), bas_at_name(:), l_sym(:) 
  integer,allocatable :: atom_num(:), block_dim(:), block_start(:) 
  complex(8) :: exp_factor
  
  pi = 4.d0 * datan(1.d0)

  ! scanner for reading user's typed comands
  ! nblocks - number of matrix blocks for the current atom in Hamiltonian matrix
  ! new_pos - position of the atom transferred to translation vector
  ! distance = distance from central atom to it's neighbour
  ! bas_at_name - name of the atom after %basis tag

!============================================================================
!
!                             Preparations
!
!============================================================================

!%%%%%%       Checking out if system.am exists       %%%%%%

  open(iuham, file = system_name, form = 'formatted', status = 'old', iostat = ierr, action = 'read' )
  if(ierr/=0) stop 'Cannot find system.am file'

!%%%%%%	        Reading system.am using tags         %%%%%%

  if( tag_exist(iuham,'cell') )then
    read(iuham,*) alat
    do i = 1,3
      read(iuham,*,iostat=ierr) ( brav_lattice(i,j), j=1,3 )
    end do
  end if 


  if ( tag_exist(iuham,'atoms') )then
    read (iuham,*) natoms
    allocate ( at_position(natoms,3), at_name(natoms) )
    do i = 1,natoms
      read (iuham,*,iostat=ierr) at_name(i), ( at_position(i,j), j=1,3 )
    end do
  end if

  if ( tag_exist( iuham , 'basis' )) then
          
          read (iuham,*) dim , nblocks
          allocate ( bas_at_name (natoms) , atom_num (natoms) , l_sym (natoms) , block_dim (natoms), block_start (natoms))
     
            do i = 1,natoms
               read (iuham,*,iostat=ierr) bas_at_name(i), atom_num(i), l_sym(i), block_dim(i) , block_start(i)
            end do
          
                !checking for duplicates 
            do i=1,natoms-1
              do j = i+1,natoms
                       if(at_name(i).eq.at_name(j))stopmark = stopmark - 2
                       if(bas_at_name(i).eq.bas_at_name(j))stopmark1 = stopmark1 -2
                       end do
            end do   
                
                if(stopmark1.ne.stopmark) stop 'Some of the basis atom name does not match with names mentioned after %atoms'
                !finished checking for duplicates

                !checking for same names       
                        do i = 1 , natoms
                          do j = 1 , natoms
                             if(at_name(i).eq.bas_at_name(j)) stopmark = stopmark + 1
                          end do
                       end do
          if(stopmark/=natoms) stop 'Some of the basis atom name does not match with names mentioned after %atoms'
                !finished checking same names
                
                !checking block_dim and block_start
            do i = 1,(natoms-1)
                if( (block_dim(i)+block_start(i)).ne.(block_start(i+1))) stop 'It must be a mistake in setting block_dim or block_start'
            end do

                !finished checking block_dim and block_start
  end if    

!       Finished reading system.am

!%%%%%%       Writing system.am data to console      %%%%%%

  write(6,*)
  write(6,*)
  write(6,*)'|system.am content|'
  write(6,*)
  write(6,*)'           bravais lattice' 
  write(6,*)
        do i = 1,3
          write(6,"(3f12.7)") (brav_lattice(i,j),j=1,3)
        enddo


  write(6,*)
  write(6,'(i2,A)') natoms,' atoms were found'
  write(6,*) 'with atomic positions'
  write(6,*)

        do i = 1, natoms
          write(6,"(i3,3(1x,f12.7))") i, ( at_position(i,j), j = 1,3 )
        enddo
 
  write(6,*)
  write(6,*)'           basis'  
  write(6,*)
  write(6,'(A,i3)')' dim:',dim 
  write(6,'(A,i3)')' nblocks:' ,nblocks
  write(6,*)  
  write(6,*) 'at_name  atom_num  l_sym   block_dim   block_start'
  do i = 1, natoms
  write(6,*) bas_at_name(i),' ', atom_num(i),' ', l_sym(i),' ', block_dim(i),' ', block_start(i)
  end do
  write(6,*)
  write(6,*)'|system.am content end|'
  write(6,*)

!       Finished writing system.am data to console

!%%%%%%              Reading Hamiltonian             %%%%%%

  write(6,*)'           Reading Hamiltonian, please wait'
  write(6,*)
  call hamiltonian_name( iuham, ham_name, ierr )
  if( ierr /= 0 )   stop 'Cannot locate hamiltonian file.'
  call read_hamiltonian( iuham, ham_name, ierr )

  select case( ierr )
    case( -14001 )
      use_hmlt = .true.
      if( rhtm == 'UNDEFINED' ) rhtm = 'TBLMTO'
    case( -14002 )
      use_hmlt = .true.
      if( rhtm == 'UNDEFINED' ) rhtm = 'ESPRESSO'
    case( :-14003, -14000:0 )
      stop 'Cannot locate or read Hamiltonian file.'
    case default
      rhtm = 'TBLMTO'
      use_hmlt = .false.
  end select

!       Finished reading Hamiltonian

!============================================================================
!
!                                 Frame
!
!============================================================================
  
  write(6,*)
  write(6,*) '  Enter number of a central atom to evaluate hopping integrals'
  write(6,*)
  read(5,*)  i_central
  
        do while(i_central>natoms.or.i_central<1) 
                write(6,*)
                        write(6,*)'   Reenter number of the central atom' 
                write(6,*)
                read(5,*) i_central
        end do

  write(6,*)
  write(6,*) '  Name and position of the central atom are'
  write(6,*)
  write(6,'(A,3f12.7)') at_name(i_central), at_position(i_central,:)

        
 
!%%%%%%      Calculating number of neighbourghs      %%%%%%

  stopmark = 0
  stopmark1 = 0

  do while(stopmark.eq.0)

  stopmark1 = 0
  nnbrs = 0

  write(6,*)
  write(6,*) '  Enter radius around the central atom in units of alat to evaluate hopping integrals'
  write(6,*)
  read(5,*)  radius

  do while( radius > 10 )
    write(6,*)
    write(6,*) ' Radius is larger than 10*alat'
    write(6,*) ' Reenter radius (or change code)'
    write(6,*)
    read(5,*) radius
  end do

  write(6,*)  

  ntransl = ceiling(radius) + 2
  
  do i = -ntransl, ntransl
    do j = -ntransl, ntransl
      do k = -ntransl, ntransl
!        translation vector is
        transl_vec = i*brav_lattice(1,:) + j*brav_lattice(2,:) + k*brav_lattice(3,:)

        do iatom = 1, natoms
          new_pos = at_position(iatom,:) + transl_vec
          distance = sqrt( ( new_pos(1) - at_position(i_central,1) )**2  +         &
                           ( new_pos(2) - at_position(i_central,2) )**2  +         &
                           ( new_pos(3) - at_position(i_central,3) )**2 )
                                
          if( distance <= radius .and. distance > 0.001 )then
            nnbrs = nnbrs + 1
            write(6,"(' Atomic positions ',3(1x,f12.7))") new_pos - at_position(i_central,:)
          end if
                                
        end do

      end do
    end do
  end do
  
  write(6,*) '  Total number of neighbours are : ', nnbrs
  write(6,*)
  write(6,*) '  Would you like to recalculate them using a different radius? [Yes/No]'
        
  do while( stopmark1.eq.0 )
    read(5,*) scanner
    scanner = adjustl(scanner)
    
    if(scanner.eq.'no'.or.scanner.eq.'No'.or.scanner.eq.'n'.or.scanner.eq.'N'.or.scanner.eq.'NO' ) then 
                        
          stopmark = 1
          stopmark1 = 1
                        
    else if( scanner .eq. 'yes' .or. scanner .eq. 'Yes' .or. scanner .eq. 'y' .or. scanner .eq. 'Y' .or. scanner .eq. 'YES' ) then
                                
          stopmark1 = 1

    else 
          write(6,*) '  Please type Yes/No'
    
    end if        
  end do   !  do while stop

  end do   !   do while
  
!     Finished calculating number of neigbours

  open( 100, file='hoppings.out', form='formatted', status='unknown', action='write' )

         
  allocate(dec_at_position(natoms,3),dec_k(3,nkp))
 
!%%%%%%	 Getting dec_at_position and dec_k matrices  %%%%%%

  dec_at_position = 0
  dec_k = 0
  brav_lattice = brav_lattice * alat
 
      do iatom = 1, natoms
          do j = 1, 3
                do k = 1, 3
                        dec_at_position(iatom,j) = at_position(iatom,k)*brav_lattice(k,j) + dec_at_position(iatom,j)
                end do
          end do
      end do


 b = 0
 M = vec_multy(brav_lattice(2,:),brav_lattice(3,:))
!     brav_lattice volume
 V = M(1)*brav_lattice(1,1) + M(2)*brav_lattice(1,2) + M(3)*brav_lattice(1,3)

      b(1,:) = M / V
      b(2,:) = vec_multy(brav_lattice(3,:),brav_lattice(1,:)) / V 
      b(3,:) = vec_multy(brav_lattice(1,:),brav_lattice(2,:)) / V 
      b = b * 2 * pi

             do ikp = 1, nkp 
                do j = 1, 3 
                   do k = 1, 3
                        dec_k(j,ikp) = dec_k(j,ikp) + bk(k,ikp)*b(k,j)
                   end do
                end do
             end do
        

!	Funished getting dec_at_position and dec_k

!%%%%%%              1 ---                                        %%%%%%
!%%%%%%   t_ij(r) = ___\     ik(R_i - R_j - transl_vec)           %%%%%%
!%%%%%%              N /   e                           * H_ij(k)  %%%%%%
!%%%%%%                ---                                        %%%%%%

  do i = -ntransl, ntransl
    do j = -ntransl, ntransl
      do k = -ntransl, ntransl

        !translation vector is
        transl_vec = i*brav_lattice(1,:) + j*brav_lattice(2,:) + k*brav_lattice(3,:)
        do iatom = 1, natoms
          new_pos = dec_at_position(iatom,:) + transl_vec
          distance = sqrt( ( new_pos(1) - dec_at_position(i_central,1) )**2  +         &
                           ( new_pos(2) - dec_at_position(i_central,2) )**2  +         &
                           ( new_pos(3) - dec_at_position(i_central,3) )**2 )
                                
          if( distance <= radius*alat .and. distance > 0.001 )then
                  
            allocate( H_r( block_dim(i_central), block_dim(iatom), nsham) )
            new_pos =  dec_at_position(i_central,:) - new_pos
          
            H_r = cmplx(0,0,8)

            do is = 1, nsham
              do ikp = 1, nkp
                k_r = dot_product( dec_k(:,ikp), -transl_vec ) 
                exp_factor = cmplx( dcos(k_r), dsin(k_r), 8 )
                             
                do m1 = block_start(i_central) + 1, block_start(i_central)+block_dim(i_central)
                  do m2 = block_start(iatom) + 1, block_start(iatom)+block_dim(iatom)
                    H_r( m1 - block_start(i_central), m2 - block_start(iatom), is ) =                                       &
                    H_r( m1 - block_start(i_central), m2 - block_start(iatom), is ) + exp_factor * h(m1-1, m2-1, ikp, is)
                  end do
                end do
              end do
            end do
            H_r = H_r / nkp

!%%% Printng out matrix of the hopping integrals to file %%%
            write(100,'(A,A,A,3f12.8)') '  Atom ', at_name(iatom), 'at position (relative to origin) ', new_pos
            write(100,'(A,f12.8,A,3f12.8)') '  distance ', distance, ' translation ', transl_vec
            write(100,*)
            do is = 1, nsham
              do m1 = 1, block_dim(i_central)
              write(100,"(30(f12.8))") ( H_r(m1,m2,is), m2 = 1, block_dim(iatom) ) 
              end do
            end do
            write(100,*)
            deallocate( H_r )
          end if
        end do
      end do
    end do
  end do


  write(6,*)
  write(6,*) 'Results have been written to hoppings.out '
  write(6,*)
  

 contains
          function vec_multy(a,b)
                            implicit none
                            real(8) :: a(3), b(3), vec_multy(3)
                                                 
                                 vec_multy(1) = a(2) * b(3) - a(3) * b(2)  
                                 vec_multy(2) = a(3) * b(1) - a(1) * b(3)  
                                 vec_multy(3) = a(1) * b(2) - a(2) * b(1)  
          end function vec_multy


 end program hk2hr
