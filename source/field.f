c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine field  --  get the potential energy functions  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "field" sets the force field potential energy functions from
c     a parameter file and modifications specified in a keyfile
c
c
      subroutine field
      use fields
      use inform
      use iounit
      use keys
      use potent
      use sizes
      implicit none
      integer i,next
      integer ia,ib
      logical header
      character*20 keyword
      character*240 record
      character*240 string

      real values(2)
      real t1, t2, t3, t4, t5, t6, t7, t8, t9
      integer nt
      real st1, st2, st3, st4, st5, st6, st7, st8, st9, st10
      data nt /0/, st1 /0/, st2 /0/, st3 /0/, st4 /0/, st5 /0/
      data st6 /0/, st7 /0/, st8 /0/, st9 /0/, st10 /0/
      save nt, st1, st2, st3, st4, st5, st6, st7, st8, st9, st10
c
c
c     set the default values for the active potentials
c
      use_bond = .true.
      use_angle = .true.
      use_strbnd = .true.
      use_urey = .true.
      use_angang = .true.
      use_opbend = .true.
      use_opdist = .true.
      use_improp = .true.
      use_imptor = .true.
      use_tors = .true.
      use_pitors = .true.
      use_strtor = .true.
      use_angtor = .true.
      use_tortor = .true.
      use_vdw = .true.
      use_repuls = .true.
      use_disp = .true.
      use_charge = .true.
      use_chgdpl = .true.
      use_dipole = .true.
      use_mpole = .true.
      use_polar = .true.
      use_chgtrn = .true.
      use_chgflx = .true.
      use_rxnfld = .false.
      use_solv = .true.
      use_metal = .false.
      use_geom = .true.
      use_extra = .true.

c     call dtime(values, t1)
c
c     read the potential energy force field parameter file
c
      call getprm
c     call dtime(values, t2)
c
c     check keywords for biopolymer atom type definitions
c
      header = .true.
      do i = 1, nkey
         next = 1
         record = keyline(i)
         call gettext (record,keyword,next)
         call upcase (keyword)
         if (keyword(1:8) .eq. 'BIOTYPE ') then
            ia = 0
            ib = 0
            string = record(next:240)
            read (string,*,err=10,end=10)  ia
            call getword (record,string,next)
            call getstring (record,string,next)
            string = record(next:240)
            read (string,*,err=10,end=10)  ib
   10       continue
            if (ia.ge.0 .and. ia.le.maxbio) then
               if (header .and. .not.silent) then
                  header = .false.
                  write (iout,20)
   20             format (/,' Additional Biopolymer Type Definitions',
     &                    //,5x,'BioType',10x,'Atom Type',/)
               end if
               biotyp(ia) = ib
               if (.not. silent) then
                  write (iout,30)  ia,ib
   30             format (1x,i8,8x,i8)
               end if

            else if (ia .gt. maxbio) then
               write (iout,40)
   40          format (/,' FIELDS  --  Too many Biopolymer Types;',
     &                    ' Increase MAXBIO')
               call fatal
            end if
         end if
      end do
c     call dtime(values, t3)
c
c     check keywords for potential function control parameters
c
      do i = 1, nkey
         record = keyline(i)
         call prmkey (record)
      end do
c     call dtime(values, t4)
      st2 = st2 + t2
      st3 = st3 + t3
      st4 = st4 + t4
      nt = nt + 1
c     write (*,*) 'field', st2/nt, st3/nt, st4/nt
      return
      end
