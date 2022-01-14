c
c
c     #########################################################
c     ##  COPYRIGHT (C) 2021 by Moses Chung & Jay W. Ponder  ##
c     ##                   All Rights Reserved               ##
c     #########################################################
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine kexind  --  exchange induction term assignment  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "kexind" assigns the size values, exponential parameter and
c     number of valence electrons for exchange induction interactions
c     and processes any new or changed values for these parameters
c
c
      subroutine kexind
      use atomid
      use atoms
      use inform
      use iounit
      use kexnd
      use keys
      use mpole
      use polpot
      use exind
      use sizes
      implicit none
      integer i,k,ii
      integer ia,ic,next
      real*8 spr,apr,epr
      logical header
      character*20 keyword
      character*240 record
      character*240 string
c
c
c     process keywords containing exchange induction parameters
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:8) .eq. 'EXCHIND ') then
            k = 0
            spr = 0.0d0
            apr = 0.0d0
            epr = 0.0d0
            call getnumb (record,k,next)
            string = record(next:240)
            read (string,*,err=10,end=10)  spr,apr,epr
   10       continue
            if (k .gt. 0) then
               if (header .and. .not.silent) then
                  header = .false.
                  write (iout,20)
   20             format (/,' Additional Exchange Induction',
     &                       ' Parameters :',
     &                    //,5x,'Atom Class',15x,'Size',11x,'Damp',
     &                       8x,'Valence'/)
               end if
               if (k .le. maxclass) then
                  peisiz(k) = spr
                  peidmp(k) = apr
                  peiele(k) = -abs(epr)
                  if (.not. silent) then
                     write (iout,30)  k,spr,apr,epr
   30                format (6x,i6,7x,2f15.4,f15.3)
                  end if
               else
                  write (iout,40)
   40             format (/,' KEXIND  --  Too many Exchange Induction',
     &                       ' Parameters')
                  abort = .true.
               end if
            end if
         end if
      end do
c
c     perform dynamic allocation of some global arrays
c
      if (allocated(sizpei))  deallocate (sizpei)
      if (allocated(dmppei))  deallocate (dmppei)
      if (allocated(elepei))  deallocate (elepei)
      allocate (sizpei(n))
      allocate (dmppei(n))
      allocate (elepei(n))
c
c     assign the repulsion size, alpha and valence parameters 
c
      nexind = n
      do i = 1, n
         sizpei(i) = 0.0d0
         dmppei(i) = 0.0d0
         elepei(i) = 0.0d0
         ic = class(i)
         if (ic .ne. 0) then
            sizpei(i) = peisiz(ic)
            dmppei(i) = peidmp(ic)
            elepei(i) = peiele(ic)
         end if
      end do
c
c     process keywords containing atom specific exchange induction
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:8) .eq. 'EXCHIND ') then
            ia = 0
            spr = 0.0d0
            apr = 0.0d0
            epr = 0.0d0
            string = record(next:240)
            read (string,*,err=70,end=70)  ia,spr,apr,epr
            if (ia.lt.0 .and. ia.ge.-n) then
               ia = -ia
               if (header .and. .not.silent) then
                  header = .false.
                  write (iout,50)
   50             format (/,' Additional Exchange Induction Values',
     &                       ' for Specific Atoms :',
     &                    //,8x,'Atom',17x,'Size',12x,'Damp',
     &                       8x,'Valence'/)
               end if
               if (.not. silent) then
                  write (iout,60)  ia,spr,apr,epr
   60             format (6x,i6,7x,2f15.4,f15.3)
               end if
               sizpei(ia) = spr
               dmppei(ia) = apr
               elepei(ia) = -abs(epr)
            end if
   70       continue
         end if
      end do
c
c     condense exchange induction sites to the list of multipole sites
c
      if (use_exind) then 
         nexind = 0
         do ii = 1, npole
            i = ipole(ii)
            if (sizpei(i) .ne. 0)  nexind = nexind + 1
            sizpei(ii) = sizpei(i)
            dmppei(ii) = dmppei(i)
            elepei(ii) = elepei(i)
         end do
      end if
c
c     turn off the exchange induction potential if not used
c
      if (nexind .eq. 0)  use_exind = .false.
      return
      end
