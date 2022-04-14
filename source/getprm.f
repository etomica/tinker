c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine getprm  --  get force field parameter file  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "getprm" finds the potential energy parameter file
c     and then opens and reads the parameters
c
c
      subroutine getprm
      use files
      use inform
      use iounit
      use keys
      use params
      implicit none
      integer i,j,iprm
      integer nask,next
      integer freeunit
      integer trimtext
      logical exist,useprm
      character*4 none
      character*20 keyword
      character*240 prmfile
      character*240 prefix
      character*240 record
      character*240 string

      real values(2)
      real t1, t2, t3, t4, t5, t6, t7, t8, t9, t10
      integer nt
      real st1, st2, st3, st4, st5, st6, st7, st8, st9, st10
      data nt /0/, st1 /0/, st2 /0/, st3 /0/, st4 /0/, st5 /0/
      data st6 /0/, st7 /0/, st8 /0/, st9 /0/, st10 /0/
      save nt, st1, st2, st3, st4, st5, st6, st7, st8, st9, st10
c
c
c     set the default name for the parameter file
c
c     call dtime(values, t1)
      useprm = .true.
      prmfile = filename(1:leng)//'.prm'
c
c     search the keyword list for the parameter filename
c
c     call dtime(values, t2)
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:11).eq.'PARAMETERS '
     &          .or. keyword(1:10).eq.'PARAMETER ') then
            string = record(next:240)
            next = 1
            call getstring (string,prmfile,next)
            if (next .eq. 1)  call gettext (string,prmfile,next)
         end if
      end do
c     call dtime(values, t3)
c
c     account for home directory abbreviation in filename
c
      if (prmfile(1:2) .eq. '~/') then
         call getenv ('HOME',prefix)
         prmfile = prefix(1:trimtext(prefix))//
     &                prmfile(2:trimtext(prmfile))
      end if
c
c     check existence of default or specified parameter file
c
      call suffix (prmfile,'prm','old')
      inquire (file=prmfile,exist=exist)
c     call dtime(values, t4)
c
c     test for user specified absence of a parameter file
c
      if (.not. exist) then
         none = prmfile(1:4)
         call upcase (none)
         if (none .eq. 'NONE') then
            exist = .true.
            useprm = .false.
         end if
      end if
c
c     try to get a parameter filename from the command line
c
      if (.not. exist) then
         call nextarg (prmfile,exist)
         if (exist) then
            call suffix (prmfile,'prm','old')
            inquire (file=prmfile,exist=exist)
         end if
      end if
c     call dtime(values, t5)
c
c     if necessary, ask for the parameter filename
c
      nask = 0
      do while (.not.exist .and. nask.lt.maxask)
         nask = nask + 1
         write (iout,10)
   10    format (/,' Enter Parameter File Name [<Enter>=NONE] :  ',$)
         read (input,20)  prmfile
   20    format (a240)
         next = 1
         call getword (prmfile,none,next)
         call upcase (none)
         if (next .eq. 1) then
            exist = .true.
            useprm = .false.
         else if (none.eq.'NONE' .and. next.eq.5) then
            exist = .true.
            useprm = .false.
         else
            if (prmfile(1:2) .eq. '~/') then
               call getenv ('HOME',prefix)
               prmfile = prefix(1:trimtext(prefix))//
     &                      prmfile(2:trimtext(prmfile))
            end if
            call suffix (prmfile,'prm','old')
            inquire (file=prmfile,exist=exist)
         end if
      end do
      if (.not. exist)  call fatal
c     call dtime(values, t6)
c
c     initialize force field control and parameter values
c
      call initprm
c
c     read the parameter file and store it for latter use
c
c     call dtime(values, t7)
      nprm = 0
      if (useprm) then
         iprm = freeunit ()
         open (unit=iprm,file=prmfile,status='old')
         rewind (unit=iprm)
         do while (.true.)
            read (iprm,30,err=50,end=50)  record
   30       format (a240)
            nprm = nprm + 1
            prmline(nprm) = record
            if (nprm .ge. maxprm) then
               write (iout,40)
   40          format (/,' GETPRM  --  Parameter File Too Large;',
     &                    ' Increase MAXPRM')
               call fatal
            end if
         end do
   50    continue
         close (unit=iprm)
      end if
c     call dtime(values, t8)
c
c     convert underbar characters to dashes in all keywords
c
      do i = 1, nprm
         next = 1
         record = prmline(i)
         call gettext (record,keyword,next)
         do j = 1, next-1
            if (record(j:j) .eq. '_')  record(j:j) = '-'
         end do
         prmline(i) = record
      end do
c     call dtime(values, t9)
c
c     get control and parameter values from the parameter file
c
      if (useprm)  call readprm
c     call dtime(values, t10)
      st2 = st2 + t2
      st3 = st3 + t3
      st4 = st4 + t4
      st5 = st5 + t5
      st6 = st6 + t6
      st7 = st7 + t7
      st8 = st8 + t8
      st9 = st9 + t9
      st10 = st10 + t10
      nt = nt + 1
c     write (*,*) 'getprm', st2/nt, st3/nt, st4/nt, st5/nt, st6/nt,
c    &                    st7/nt, st8/nt, st9/nt, st10/nt
      return
      end
