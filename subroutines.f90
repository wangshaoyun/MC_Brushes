module subroutines
use global_variables
use compute_energy
implicit none
contains

!################Initialize Position################!
!origin is on a corner of square
!the former 3 columns is positions, and the 4th column is charge
!the charge of end particle is qq
!the free particle is -qq

subroutine Initialize(pos)
use global_variables
implicit none

!!!!!!!!data dictionary!!!!!!!!
	real*8, dimension(NN*Nz,4), intent(out):: pos
  real*8 :: dl                		 !distance between two particle
  integer:: i,j,k,l,n                        !No.l particle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	pos=0.
  
!链上粒子初始化，立方晶格
	n=0
	l=NN**(1./2)      ! l: number of chains in each direction
  dl=LL/l         ! dl: distance between two chains
  do i=1,l
    do j=1,l
      do k=1,Nz-Nq
        n=n+1
        pos(n,1)=(i-0.5)*dl
        pos(n,2)=(j-0.5)*dl
        pos(n,3)=(k-1)*dz
        if (mod(n, Nz-Nq)==0) then
          pos(n,4)=qq
        end if
      end do
    end do
  end do
!游离粒子初始化，放到链末端上方两个化学键长处
  do i=1,l
    do j=1,l
      do k=1, Nq
        n=n+1
        pos(n,1)=(i-0.01)*dl
        pos(n,2)=(j-0.01)*dl
        pos(n,3)=(k+1)*dz!+Lz/2
        pos(n,4)=-qq/abs(qq)
      end do
    end do
  end do

end subroutine Initialize
!##############End Initialize Position##############!


subroutine Initialize_EC_vector(pos, EC_vector)
use global_variables
implicit none
!!!!!!!!data dictionary!!!!!!!!
	real*8, dimension(NN*Nz,4), intent(in):: pos
  real*8, dimension(Kfourier+3,5), intent(out) :: EC_vector
  real*8 :: rr                		 !distance between two particle
  integer:: i,j,k,n,p,q,s,t                        !No.l particle
  real*8, dimension(3) :: rpq
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  EC_vector=0.
!计算库伦势的四部分

  do s=1,NN*(Nq+1)    !self energy
    if(s<=NN) then
      p=(Nz-Nq)*s
    else 
      p=NN*(Nz-Nq)+(s-NN)
    end if
    EC_vector(1,1)=EC_vector(1,1)+alpha/sqrt(pi)*pos(p,4)**2
  end do

  do s=1,NN*(Nq+1)    !slab energy
    if(s<=NN) then
      p=(Nz-Nq)*s
    else 
      p=NN*(Nz-Nq)+(s-NN)
    end if
    EC_vector(2,1)=EC_vector(2,1)+pos(p,4)*pos(p,3)
  end do

  do s=1,NN*(Nq+1)  !(Nq+1) is the number of charge
    do t=1,NN*(Nq+1)
      if(s==t) cycle
!transfer s,t to p,q
      if(s<=NN) then
        p=(Nz-Nq)*s
      else 
        p=NN*(Nz-Nq)+(s-NN)
      end if
      if(t<=NN) then
        q=(Nz-Nq)*t
      else 
        q=NN*(Nz-Nq)+(t-NN)
      end if
      rpq(1:3)=pos(q,1:3)-pos(p,1:3)
      rr=sqrt(dot_product(rpq, rpq))
      EC_vector(3,1)=EC_vector(3,1)+pos(p,4)*pos(q,4)*erfc(alpha*rr)/rr/2
    end do
  end do

  n=0
  do i=-Kmax, Kmax
    do j=-Kmax, Kmax
      do k=-Kmax, Kmax
        if((i*i+j*j+k*k)>(Kmax*Kmax) .or. (i*i+j*j+k*k)==0) cycle
        n=n+1
        EC_vector(n+3,1)=2*pi*i/LL
        EC_vector(n+3,2)=2*pi*j/LL
        EC_vector(n+3,3)=2*pi*k/(Lz*6)  !empty space incerted is 5Lz
        do s=1, NN*(Nq+1)
          if(s<=NN) then
            p=(Nz-Nq)*s
          else 
            p=NN*(Nz-Nq)+(s-NN)
          end if
          EC_vector(n+3,4)=EC_vector(n+3,4)+pos(p,4)               &
             *cos(dot_product(EC_vector(n+3,1:3),pos(p,1:3)))
          EC_vector(n+3,5)=EC_vector(n+3,5)+pos(p,4)               &
             *sin(dot_product(EC_vector(n+3,1:3),pos(p,1:3)))
        end do
      end do
    end do
  end do
!    write(*,*) n, Kfourier
end subroutine Initialize_EC_vector


!####################verlet list####################!

subroutine verlet_lj_list(lj_list, pos)
use global_variables
implicit none
integer, dimension(NN*Nz+1, NN*Nz), intent(out):: lj_list
real*8, dimension(NN*Nz,4), intent(in) :: pos
integer :: i,j
real*8, dimension(3) :: rij
real*8 :: rr
lj_list=0

do i=1, NN*Nz-1
  do j=i+1, NN*Nz
    rij=pos(i,1:3)-pos(j,1:3)
    rij(1)=((rij(1)+LL/2)-LL*floor((rij(1)+LL/2)/LL))-LL/2
		rij(2)=((rij(2)+LL/2)-LL*floor((rij(2)+LL/2)/LL))-LL/2
    rr=dot_product(rij, rij)
    if (rr<(rv_lj*rv_lj)) then
       lj_list(NN*Nz+1,i)=lj_list(NN*Nz+1,i)+1
       lj_list(NN*Nz+1,j)=lj_list(NN*Nz+1,j)+1
       lj_list(lj_list(NN*Nz+1,i), i)=j
       lj_list(lj_list(NN*Nz+1,j), j)=i
    end if
  end do
end do
end subroutine verlet_lj_list
!##################End verlet list##################!


!####################verlet list####################!

subroutine verlet_real_list(real_list, pos)
use global_variables
implicit none
integer, dimension(NN*(Nq+1)+1, NN*(Nq+1)), intent(out) :: real_list
real*8, dimension(NN*Nz, 4), intent(in) :: pos
integer :: i, j, p, q, r
real*8, dimension(3) :: rij
real*8 :: rr

real_list=0

  do p=1, NN*(Nq+1)-1
    do q=p+1, NN*(Nq+1)
      if(p<=NN) then
        i=(Nz-Nq)*p
      else
        i=NN*(Nz-Nq)+(p-NN)
      end if
      if(q<=NN) then
        j=(Nz-Nq)*q
      else
        j=NN*(Nz-Nq)+(q-NN)
      end if
      rij=pos(i,1:3)-pos(j,1:3)
      rij(1)=((rij(1)+LL/2)-LL*floor((rij(1)+LL/2)/LL))-LL/2
		  rij(2)=((rij(2)+LL/2)-LL*floor((rij(2)+LL/2)/LL))-LL/2
      rr=dot_product(rij, rij)
      if (rr<(rv_real*rv_real)) then
        real_list(NN*(Nq+1)+1,p)=real_list(NN*(Nq+1)+1,p)+1
        real_list(NN*(Nq+1)+1,q)=real_list(NN*(Nq+1)+1,q)+1
        real_list(real_list(NN*(Nq+1)+1, p), p)=j
        real_list(real_list(NN*(Nq+1)+1, q), q)=i
      end if
    end do
  end do

end subroutine verlet_real_list
!##################End verlet list##################!


!##########Choose particle which will move##########!
!不能选择底端固定的粒子

subroutine Choose_Particle(l, pos)
use global_variables
implicit none
!!!!!!!!data dictionary!!!!!!!!
	integer, intent(out):: l          !choosed particle
	real*8, dimension(NN*Nz,4), intent(in) :: pos
	real*8 :: rnd                       !random number 0-1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	l=1
!不能选择底端粒子移动
  do while((mod(l,(Nz-Nq))==1 .and. l<=NN*(Nz-Nq)))
		call random_number(rnd)
		l=ceiling(rnd*NN*Nz)
	end do

end subroutine Choose_Particle
!########End Choose particle which will move########!


!####################New Position###################!
!移动被选中粒子，产生三个随机数得到三方向位移，得到新位置分布，注意z方向没有周期性边界条件

subroutine New_Position(pos, pos1, EC_vector, EC_vector1, real_list, l)
use global_variables
implicit none
!!!!!!!!data dictionary!!!!!!!!
	real*8, dimension(NN*Nz,4), intent(in):: pos
	real*8, dimension(NN*Nz,4), intent(out):: pos1
  real*8, dimension(Kfourier+3,5), intent(in) :: EC_vector
  real*8, dimension(Kfourier+3,5), intent(out) :: EC_vector1
  integer, dimension(NN*(Nq+1)+1, NN*(Nq+1)), intent(in):: real_list
	integer, intent(in):: l
	integer :: n,p,s
	real*8, dimension(3) :: rlp, rlp1
	real*8 :: rndx, rndy, rndz, rr, rr1, a, b, k2, rhok2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  pos1=pos
  EC_vector1=EC_vector
  
  call random_number(rndx)
	call random_number(rndy)
	call random_number(rndz)
	pos1(l,1)=pos(l,1)+(rndx-0.5)*dr
	pos1(l,2)=pos(l,2)+(rndy-0.5)*dr
	pos1(l,3)=pos(l,3)+(rndz-0.5)*dr
  pos1(l,1)=pos1(l,1)-LL*floor(pos1(l,1)/LL)
	pos1(l,2)=pos1(l,2)-LL*floor(pos1(l,2)/LL)


if( pos(l,4) /=0 ) then

  EC_vector1(2,1)=EC_vector(2,1)-pos(l,3)*pos(l,4)+pos1(l,3)*pos1(l,4)

  do s=1, NN*(Nq+1)
    if (s<=NN) then
      p=s*(Nz-Nq);
    else 
      p=NN*(Nz-Nq)+(s-NN)
    end if
    if( l==p ) cycle
    rlp(1:3)=pos(p,1:3)-pos(l,1:3)
    rlp1(1:3)=pos1(p,1:3)-pos1(l,1:3)
    rr=sqrt(dot_product(rlp, rlp))
    rr1=sqrt(dot_product(rlp1, rlp1))
    EC_vector1(3,1)=EC_vector1(3,1)-pos(l,4)*pos(p,4)*erfc(alpha*rr)/rr/2      &
&               +pos1(l,4)*pos1(p,4)*erfc(alpha*rr1)/rr1/2      
  end do

  do n=4, Kfourier+3
    EC_vector1(n,4)=EC_vector1(n,4)-pos(l,4)                                   &
          *cos(dot_product(EC_vector(n,1:3), pos(l,1:3)))                      &
          +pos1(l,4)*cos(dot_product(EC_vector(n,1:3),pos1(l,1:3)))
    EC_vector1(n,5)=EC_vector1(n,5)-pos(l,4)                                   &
          *sin(dot_product(EC_vector(n,1:3),pos(l,1:3)))                       &
          +pos1(l,4)*sin(dot_product(EC_vector(n,1:3),pos1(l,1:3)))
  end do

end if

end subroutine New_Position
!##################End New Position#################!


!##################Move particle or not#################!
!能量差小于0，则移动，若大于0，则产生随机数，小于玻尔兹曼因子则移动

subroutine Move_or_not(pos, pos1, EC_vector, EC_vector1, l, DeltaE)
use global_variables
implicit none
!!!!!!!!data dictionary!!!!!!!!
	real*8, dimension(NN*Nz,4), intent(inout):: pos
	real*8, dimension(NN*Nz,4), intent(in):: pos1
  real*8, dimension(Kfourier+3,5), intent(inout) :: EC_vector
  real*8, dimension(Kfourier+3,5), intent(in) :: EC_vector1
	real*8, intent(in) :: DeltaE
  integer, intent(in):: l
	real*8  :: rnd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	if (DeltaE<0) then
		pos=pos1
		EC_vector=EC_vector1
	else 
		call random_number(rnd)
		if (rnd<Exp(-Beta*DeltaE)) then
			pos=pos1
			EC_vector=EC_vector1
		endif
	endif

end subroutine Move_or_not
!################End Move particle or not################!


!######Height Distribution funciton of End particle#######!

subroutine Height_Dist_Function(pos, HistHdfEnd1, HistHdfEnd2, HistHdfEnd3)
!!!!!!!!数据词典!!!!!!!!
	real*8, dimension(NN*Nz,4), intent(in) :: pos
	real*8, dimension(SizeHist), intent(inout) :: HistHdfEnd1
  real*8, dimension(SizeHist), intent(inout) :: HistHdfEnd2
  real*8, dimension(SizeHist), intent(inout) :: HistHdfEnd3
  real*8 :: rz
	integer :: i,j,k
!!!!!!!!!!!!!!!!!!!!!!!

! hdf of charged particle in the end
  do i=1,NN
			j=i*(Nz-Nq)
	    rz=pos(j,3)
      k=ceiling(rz/deltaR)
      if (k<=0 .or. k>=SizeHist) then
        write(*,*) 'hdf of end charged particle, k==0 or k>SizeHist',k,rz,j
        cycle
      end if
      HistHdfEnd1(k)=HistHdfEnd1(k)+1
  end do

!hdf of chains
  do i=1,NN*(Nz-Nq)
      if (mod(i,(Nz-Nq))==1) cycle
	    rz=pos(i,3)
      k=ceiling(rz/deltaR)
      if (k<=0 .or. k>=SizeHist) then
        write(*,*) 'hdf of chains, k==0 or k>SizeHist',k
        cycle
      end if
      HistHdfEnd2(k)=HistHdfEnd2(k)+1
  end do

!hdf of counter ions
  do i=NN*(Nz-Nq)+1, NN*Nz
	    rz=pos(i,3)
      k=ceiling(rz/deltaR)
      if (k<=0 .or. k>=SizeHist) then
        write(*,*) 'hdf of counter ions, k==0 or k>SizeHist',k
        cycle
      end if
      HistHdfEnd3(k)=HistHdfEnd3(k)+1
  end do


end subroutine Height_Dist_Function
!####End Height Distribution funciton of End particle#####!


!#################End to Center distance##################!

subroutine end_center_distance(pos, RG, RGz, R2, Rz2, height)
use global_variables
implicit none
real*8, dimension(NN*Nz,4), intent(in) :: pos
real*8, intent(out) :: RG
real*8, intent(out) :: RGz
real*8, intent(out) :: R2
real*8, intent(out) :: Rz2
real*8, intent(out) :: height
real*8, dimension(3) :: RR
integer :: i,j,k

 RG=0.
 RGz=0.
 R2=0.
 Rz2=0.
 height=0.

do i=1,NN

  do j=1, Nz-Nq
    do k=1, Nz-Nq
      if (j==k) cycle
      RR=pos((i-1)*(Nz-Nq)+j,1:3)-pos((i-1)*(Nz-Nq)+k,1:3)
      RR(1)=((RR(1)+LL/2)-LL*floor((RR(1)+LL/2)/LL))-LL/2
		  RR(2)=((RR(2)+LL/2)-LL*floor((RR(2)+LL/2)/LL))-LL/2
      RG=RG+dot_product(RR, RR)/NN/Nz/Nz/2
      RGz=RGz+RR(3)*RR(3)/NN/Nz/Nz/2
    end do
  end do
  
  RR=pos(i*(Nz-Nq),1:3)-pos((i-1)*(Nz-Nq)+1,1:3)
  RR(1)=((RR(1)+LL/2)-LL*floor((RR(1)+LL/2)/LL))-LL/2
  RR(2)=((RR(2)+LL/2)-LL*floor((RR(2)+LL/2)/LL))-LL/2
  R2=R2+dot_product(RR,RR)/NN
  Rz2=Rz2+RR(3)*RR(3)/NN
  
  do j=1,Nz
    height=height+pos((i-1)*(Nz-Nq)+j,3)/NN/Nz
  end do
  
end do

end subroutine end_center_distance
!#################End to Center distance##################!


!################Output Physical Quantity#################!

subroutine write_pos0(pos)
implicit none
	real*8, dimension(NN*Nz,4), intent(in) :: pos
	open(10, file='/home/a/wangshaoyun/data/20171204/2/pos0.txt')
    write(10,110) transpose(pos)
    110 format(4F15.6)
  close(10)
end subroutine write_pos0

subroutine write_pos1(pos, steps)
implicit none
	real*8, dimension(NN*Nz,4), intent(in) :: pos
	integer, intent(in) ::steps
	open(20, file='/home/a/wangshaoyun/data/20171204/2/pos1.txt')
	  write(20,*) steps
    write(20,120) transpose(pos)
    120 format(4F15.6)
  close(20)
end subroutine write_pos1

subroutine write_energy(EE, DeltaE)
implicit none
	real*8, intent(in) :: EE
  real*8, intent(in) :: DeltaE
  real*8 :: pb
  pb=1.
  if (DeltaE>0) then
    pb=exp(-Beta*DeltaE)
  end if
	open(UNIT=30, Position='Append',file='/home/a/wangshaoyun/data/20171204/2/energy.txt')
    write(30,130) EE, DeltaE, pb
    130 format(3F15.6)
  close(30)
end subroutine write_energy

subroutine write_hdf(HistHdfEnd1, HistHdfEnd2, HistHdfEnd3)
use global_variables
implicit none
	real*8, dimension(SizeHist), intent(out) :: HistHdfEnd1
	real*8, dimension(SizeHist), intent(out) :: HistHdfEnd2
	real*8, dimension(SizeHist), intent(out) :: HistHdfEnd3
	integer :: i
	open(14, file='/home/a/wangshaoyun/data/20171204/2/hdf.txt')
  do i=1, SizeHist
    write(14,140) i*Lz/SizeHist, HistHdfEnd1(i), HistHdfEnd2(i), HistHdfEnd3(i)
    140 format(4F15.6)
  end do
  close(14)
end subroutine write_hdf

subroutine write_rec(RG, RGz, R2, Rz2, height)
use global_variables
implicit none
  real*8, intent(inout) :: RG
  real*8, intent(inout) :: RGz
  real*8, intent(inout) :: R2
  real*8, intent(inout) :: Rz2
  real*8, intent(inout) :: height
  open(UNIT=15, Position='Append', file='/home/a/wangshaoyun/data/20171204/2/rec.txt')
    write(15,150) Rho, Nz*1., RG, RGz, R2, Rz2, height
    150 format(7F15.6)
  close(15)
end subroutine write_rec

!###############End Output Physical Quantity##############!


end module subroutines




































