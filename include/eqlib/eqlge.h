! eqlge.h
!     Data specific to the generic ion-exchange model.
character(len=56) :: ugexr
character(len=24) :: ugexmo
character(len=24) :: ugexp
character(len=24) :: ugexs
character(len=24) :: ugexsr
character(len=8) :: ugexj
character(len=8) :: uhfgex
character(len=8) :: uvfgex
character(len=8) :: uxkgex

common /eqlgec/ ugexj(jet_par,net_par),ugexmo(net_par),ugexp(net_par),ugexr(iet_par,jet_par,net_par),ugexs(iet_par,jet_par,net_par),ugexsr(2,iet_par,jet_par,net_par),uhfgex(iet_par,jet_par,net_par),uvfgex(iet_par,jet_par,net_par),uxkgex(iet_par,jet_par,net_par)

integer :: jern1
integer :: jern2
integer :: jgext
integer :: kern1
integer :: kern2
integer :: kgexsa
integer :: ngexro
integer :: ngexrt
integer :: ngexsa
integer :: ngexso
integer :: ngext

common /eqlgei/ jern1(jet_par,net_par),jern2(jet_par,net_par),jgext(net_par),kern1(net_par),kern2(net_par),kgexsa(ket_par,net_par),ngexro(iet_par,jet_par,net_par),ngexrt(jet_par,net_par),ngexsa(iet_par,jet_par,net_par),ngexso(iet_par,jet_par,net_par),ngext(jet_par,net_par)

real(kind=8) :: cgexj
real(kind=8) :: mgext
real(kind=8) :: mwtges
real(kind=8) :: tgexp
real(kind=8) :: xhfgex
real(kind=8) :: xlkgex
real(kind=8) :: xvfgex
real(kind=8) :: zgexj

common /eqlge/ cgexj(jet_par,net_par),mgext(jet_par,net_par),mwtges(net_par),tgexp(net_par),xhfgex(iet_par,jet_par,net_par),xlkgex(iet_par,jet_par,net_par),xvfgex(iet_par,jet_par,net_par),zgexj(jet_par,net_par)

