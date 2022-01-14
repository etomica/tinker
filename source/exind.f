c
c
c     #########################################################
c     ##  COPYRIGHT (C) 2021 by Moses Chung & Jay W. Ponder  ##
c     ##                   All Rights Reserved               ##
c     #########################################################
c
c     ##################################################################
c     ##                                                              ##
c     ##  module exind  --  exchange induction for current structure  ##        
c     ##                                                              ##
c     ##################################################################
c
c
c     nexind    total number of exchange induction sites in the system
c     sizpei    Pauli repulsion size parameter value at each site
c     dmppei    Pauli repulsion alpha damping value at each site
c     elepei    Pauli repulsion valence electrons at each site
c
c
      module exind
      implicit none
      integer nexind
      real*8, allocatable :: sizpei(:)
      real*8, allocatable :: dmppei(:)
      real*8, allocatable :: elepei(:)
      save
      end
