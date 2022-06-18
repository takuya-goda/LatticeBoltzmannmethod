program lbm
   !移動円柱流れ
   use :: iso_fortran_env
   implicit none

   !パラメーターの設定
   integer, parameter :: NumCell_x = 500
   integer, parameter :: NumCell_y = NumCell_x
   integer, parameter :: Nx = NumCell_x
   integer, parameter :: Ny = NumCell_y

   
   !角柱の設定
   integer, parameter :: prizm_top    = NumCell_x/2 + 10
   integer, parameter :: prizm_bottom = NumCell_x/2 - 10
   integer, parameter :: prizm_left   = NumCell_y/2 + 10
   integer, parameter :: prizm_right  = NumCell_y/2 - 10

   real(real64), parameter :: rho0 = 1.0d0
   real(real64), parameter :: rho1 = 1.0d0
   real(real64), parameter :: u0 = 0.03d0
   real(real64), parameter :: h0 = dble(Nx)
   real(real64), parameter :: Re = 1000d0
   real(real64), parameter :: nu = u0*h0/Re
   real(real64), parameter :: dp = 12d0*rho0*nu/h0**2*u0

   real(real64), parameter :: dx = 1d0/dble(NumCell_x)
   real(real64), parameter :: ds = 1d0
   real(real64), parameter :: Uwall = 0.1d0
   real(real64), parameter :: dt = Uwall*dx
   real(real64), parameter :: kineticViscosity = Uwall/dx/Re !動粘性係数

   real(real64), parameter :: RelaxTime = 3d0*kineticViscosity + 0.5d0 !(1.49)式

   ! LBMの計算で使用するパラメーター設定
   integer, parameter :: Center = 0
   integer, parameter :: Right = 1
   integer, parameter :: Up = 2
   integer, parameter :: Left = 3
   integer, parameter :: Down = 4
   integer, parameter :: UpRight = 5
   integer, parameter :: UpLeft = 6
   integer, parameter :: DownLeft = 7
   integer, parameter :: DownRight = 8
   integer, parameter :: First = Center
   integer, parameter :: Last = DownRight

   ! すべり無し
   integer, parameter :: Opposite(First:Last) = (/Center, Left, Down, Right, Up, &
                                                  DownLeft, DownRight, UpRight, Upleft/)

   ! 完全すべり
   integer, parameter :: OppositeIB(First:Last) = (/Center, Left, Down, Right, Up, &
                                                    DownRight, Downleft, UpLeft, UpRight/)

   real(real64), parameter :: Weight(First:Last) = (/4d0/9d0, &
                                                1d0/9d0, 1d0/9d0, 1d0/9d0, 1d0/9d0, &
                                                1d0/36d0, 1d0/36d0, 1d0/36d0, 1d0/36d0/)

   integer, parameter :: ConvVelx(First:Last) = (/0, 1, 0, -1, 0, 1, -1, -1, 1/)
   integer, parameter :: ConvVely(First:Last) = (/0, 0, 1, 0, -1, 1, 1, -1, -1/)

   !柱の設定
   real(real64), parameter :: xc = 0.5d0*dble(Nx) ! 円柱中心
   real(real64), parameter :: yc = 0.5d0*dble(Ny)
   real(real64), parameter :: D = 50d0            ! 円柱直径
   real(real64), parameter :: epsilon = 1.0e-14
   real(real64), parameter :: pi = dacos(-1.d0)	!円周率
   integer, parameter :: nb = 4*(int(D) + 1) ! 境界点数
   real(real64), parameter :: dv = pi*D/(dble(nb))*ds

   real(real64), dimension(NumCell_x, NumCell_x) ::velx
   real(real64), dimension(NumCell_x, NumCell_x) ::vely
   real(real64), dimension(NumCell_x, NumCell_x) ::dens
   real(real64), dimension(First:Last, NumCell_x, NumCell_x) ::f
   real(real64), dimension(First:Last, NumCell_x, NumCell_x) ::f_eq

   real(real64) :: xc_move,yc_move

   integer, parameter :: Nt = 100000
   integer, parameter :: outputstep = 500
   integer, parameter :: itr_max = 5

   integer :: n, direction, n2
   integer :: i, j, l
   integer :: i_test,j_test

   i_test = 254
   j_test = 223

   ! 初期値設定
   f(:, :, :) = 0d0
   velx(:, :) = 0d0
   vely(:, :) = 0d0
   dens(:, :) = rho0
   xc_move = xc
   yc_move = yc

   print *, RelaxTime   

   do n = 1, Nt
      ! fの計算
      call computeEquilibrium(f_eq, velx, vely, dens)
      call collide(f, f_eq)
      call stream(f)
      call imposeBoundary(f)
      call immersedBoundaryMethod(f,xc_move,yc_move,n) !円柱
      call computef(f, velx, vely, dens)      

      ! 出力
      if (mod(n, outputstep) == 0) then
         call output(dens,velx,vely,n,xc_move,yc_move)
      ! elseif (mod(n, 1) == 0 .and. n < 500) then
      !    call output(dens,velx,vely,n)
      end if
   end do

contains

   subroutine output(dens, velx, vely, time, xc_move, yc_move)
      implicit none
      real(real64), intent(in) :: dens(1:nx, 1:ny)
      real(real64), intent(in) :: velx(1:nx, 1:ny)
      real(real64), intent(in) :: vely(1:nx, 1:ny)
      integer, intent(in) :: time
      real(real64), intent(in) :: xc_move
      real(real64), intent(in) :: yc_move

      character(128) filename1
      character(128) filename2
      integer :: i, j

      write (filename1, '("9.output/dens",i5.5,".txt")') time

      open (unit=10, file=filename1)
      do j = 1, ny
      do i = 1, nx
         if((i-xc_move)**2 + (j-yc_move)**2 <=(0.5d0*D)**2) then
            write (10, '(2i8,3f23.14)') i, j, 1.5d0, 0d0, 0d0
         else
            write (10, '(2i8,3f23.14)') i, j, dens(i, j), velx(i, j), vely(i, j)
         endif
      end do
      end do
      close (10)

   end subroutine

   subroutine computeEquilibrium(f_eq, velx, vely, dens)
      implicit none
      real(real64), intent(inout)  :: f_eq(First:Last, 1:Nx, 1:Ny)
      real(real64), intent(in)     :: velx(1:Nx, 1:Ny)
      real(real64), intent(in)     :: vely(1:Nx, 1:Ny)
      real(real64), intent(in)     :: dens(1:Nx, 1:Ny)

      real(real64) ::u, v, conv_velo, velo_square, pressure
      integer ::i, j, direction

      do j = 1, Ny
      do i = 1, Nx
         u = velx(i, j)
         v = vely(i, j)
         ! if (i==1) u = uwall !境界条件
         velo_square = u*u + v*v
         do direction = First, Last
            conv_velo = u*Convvelx(direction) + v*Convvely(direction)
            f_eq(direction, i, j) = weight(direction)*dens(i, j) &
                                    *(1d0 + 3d0*conv_velo + 4.5d0*conv_velo*conv_velo - 1.5d0*velo_square)
         end do
      end do
      end do

   end subroutine

   subroutine collide(f, f_eq)
      implicit none
      real(real64), intent(inout) :: f(First:Last, 1:Nx, 1:Ny)
      real(real64), intent(in)    :: f_eq(First:Last, 1:Nx, 1:Ny)
      integer :: i, j

      do j = 1, Ny
      do i = 1, Nx
         f(:, i, j) = f(:, i, j) + (f_eq(:, i, j) - f(:, i, j))/RelaxTime
      end do
      end do

   end subroutine

   subroutine stream(f)
      implicit none
      real(real64), intent(inout)::f(First:Last, 1:Nx, 1:Ny)
      integer ::i, j

      do j = 1, Ny
         do i = Nx, 2, -1 !RIGHT TO LEFT​
            f(Right, i, j) = f(Right, i - 1, j)
         end do
         do i = 1, Nx - 1 !LEFT TO RIGHT​
            f(Left, i, j) = f(Left, i + 1, j)
         end do
      end do

      do j = Ny, 2, -1 !TOP TO BOTTOM​
         do i = 1, Nx
            f(Up, i, j) = f(Up, i, j - 1)
         end do
         do i = Nx, 2, -1
            f(UpRight, i, j) = f(UpRight, i - 1, j - 1)
         end do
         do i = 1, Nx - 1
            f(UpLeft, i, j) = f(UpLeft, i + 1, j - 1)
         end do
      end do

      do j = 1, Ny - 1 !BOTTOM TO TOP​
         do i = 1, Nx
            f(Down, i, j) = f(Down, i, j + 1)
         end do
         do i = 1, Nx - 1
            f(DownLeft, i, j) = f(DownLeft, i + 1, j + 1)
         end do
         do i = Nx, 2, -1
            f(DownRight, i, j) = f(DownRight, i - 1, j + 1)
         end do
      end do

   end subroutine

   subroutine imposeBoundary(f)
      implicit none
      real(real64), intent(inout)::f(First:Last, 1:Nx, 1:Ny)
      integer :: i, j

      do j = 1, Ny
         !bounce back on west boundary​
         f(Right     , 1, j) = f(Opposite(Right    ), 1, j)
         f(UpRight   , 1, j) = f(Opposite(UpRight  ), 1, j)
         f(DownRight , 1, j) = f(Opposite(DownRight), 1, j)
         !bounce back on east boundary​
         f(Left     , Nx, j) = f(Opposite(Left    ), Nx, j)
         f(DownLeft , Nx, j) = f(Opposite(DownLeft), Nx, j)
         f(UpLeft   , Nx, j) = f(Opposite(UpLeft  ), Nx, j)
      end do

      !bounce back on south boundary​
      do i = 1, Nx
         f(Up     , i, 1) = f(OppositeIB(Up     ), i, 1)
         f(UpRight, i, 1) = f(OppositeIB(UpRight), i, 1)
         f(UpLeft , i, 1) = f(OppositeIB(UpLeft ), i, 1)
      end do

      !bounce back on nouth boundary​
      do i = 1, Nx
         f(Down     , i, Ny) = f(OppositeIB(Down     ), i, Ny)
         f(DownRight, i, Ny) = f(OppositeIB(DownRight), i, Ny)
         f(DownLeft , i, Ny) = f(OppositeIB(DownLeft ), i, Ny)
      end do


   end subroutine

   subroutine immersedBoundaryMethod(f,xc_move,yc_move,n)
      implicit none
      real(real64), intent(inout)::f(First:Last, 1:Nx, 1:Ny)
      real(real64),intent(inout) :: xc_move
      real(real64),intent(inout) :: yc_move
      integer,intent(in) :: n

      real(real64) :: xb(1:nb)
      real(real64) :: yb(1:nb)
      real(real64) :: ub(1:nb)
      real(real64) :: vb(1:nb)
      real(real64) :: ud(1:nb)
      real(real64) :: vd(1:nb)
      real(real64) :: gbx(1:nb)
      real(real64) :: gby(1:nb)

      real(real64) :: us(1:nx,1:ny)
      real(real64) :: vs(1:nx,1:ny)
      real(real64) :: gx(1:nx,1:ny)
      real(real64) :: gy(1:nx,1:ny)
      real(real64) :: uitr(1:nx,1:ny)
      real(real64) :: vitr(1:nx,1:ny)

      integer i,j,k,l,itr
      real(real64) xc_velo,xc_old
      integer nxmin,nxmax,nymin,nymax

      ! 円柱の移動
      xc_old = xc_move
      xc_move = xc + 100d0*dcos(pi/2d0 + pi*dble(n)/50000d0)
      xc_velo = (xc_move- xc_old)/RelaxTime
      print *, n,xc + 10d0*dcos(pi/2d0 + pi*dble(n)/50000d0), &
               (xc_move- xc_old)/RelaxTime
      nxmin = floor(xc - D*0.5d0 - 2.0d0)
      nxmax = ceiling(xc + D*0.5d0 + 2.0d0)
      nymin = floor(yc - D*0.5d0 - 2.0d0)
      nymax = ceiling(yc + D*0.5d0 + 2.0d0)

      ! ===境界点の定義
      do k=1,nb
         xb(k) = xc_move + 0.5d0*D*dcos(2d0*pi/dble(nb)*k) !cos(2/pi)=>90度
         yb(k) = yc_move + 0.5d0*D*dsin(2d0*pi/dble(nb)*k)


         ! 境界で満たすべき流速 P57最後、3.6式で使用
         ! ud(k) = 0d0
         ud(k) = xc_velo
         vd(k) = 0d0
      enddo

      ! ===固体付近のセルの流速の計算(仮)
      do j=nymin,nymax
      do i=nxmin,nxmax
         us(i,j) = f(Right,i,j) - f(Left,i,j) + f(UpRight,i,j) - f(UpLeft,i,j)    + f(DownRight,i,j) - f(DownLeft,i,j)
         vs(i,j) = f(Up,i,j)    - f(Down,i,j) + f(UpRight,i,j) - f(DownRight,i,j) + f(UpLeft,i,j)    - f(DownLeft,i,j)
      enddo
      enddo

      ! ===境界点上の流速の計算（仮）
      do k=1,nb
         ub(k) = 0d0
         vb(k) = 0d0
         !　境界点上の±２の範囲に入るセルが対象
         ! 3.3
         do j=floor(yb(k)-2d0),ceiling(yb(k)+2d0)
         do i=floor(xb(k)-2d0),ceiling(xb(k)+2d0)
            ub(k) = ub(k) + us(i,j)*w(i-xb(k))*w(j-yb(k))*ds**2
            vb(k) = vb(k) + vs(i,j)*w(i-xb(k))*w(j-yb(k))*ds**2
         enddo
         enddo
      enddo

      ! ===体積力の計算
      ! step0 境界点上の体積力の算出 3.6
      do k=1,nb
         gbx(k) = (ud(k)-ub(k))/ds
         gby(k) = (vd(k)-vb(k))/ds
      enddo

      ! ===反復計算
      do itr = 1,itr_max
         ! step1 格子点上の体積力の計算
         gx(:,:) = 0d0
         gy(:,:) = 0d0
         ! 3.7式
         do k=1,nb
            do j=floor(yb(k)-2d0),ceiling(yb(k)-2d0)
            do i=floor(xb(k)-2d0),ceiling(xb(k)-2d0)
               gx(i,j) = gx(i,j) + gbx(k)*w(i-xb(k))*w(j-yb(k))*dv
               gy(i,j) = gy(i,j) + gby(k)*w(i-xb(k))*w(j-yb(k))*dv
            enddo
            enddo
         enddo

         ! step2 格子点上の流速の修正
         do j=nymin,nymax
         do i=nxmin,nxmax
            uitr(i,j) = us(i,j) + ds*gx(i,j)
            vitr(i,j) = vs(i,j) + ds*gy(i,j)
         enddo
         enddo

         ! step3 境界点上の流速の内挿
         do k=1,nb
            ub(k) = 0d0
            vb(k) = 0d0
            do j=floor(yb(k)-2d0),ceiling(yb(k)-2d0)
            do i=floor(xb(k)-2d0),ceiling(xb(k)-2d0)
               ub(k) = ub(k) + uitr(i,j)*w(i-xb(k))*w(j-yb(k))*ds**2
               vb(k) = vb(k) + vitr(i,j)*w(i-xb(k))*w(j-yb(k))*ds**2
            enddo
            enddo
         enddo

         ! step4 境界点上ですべり無し条件を満たすように体積力を修正
         ! 3.10
         do k=1,nb
            gbx(k) = gbx(k) + (ud(k)-ub(k))/ds
            gby(k) = gby(k) + (vd(k)-vb(k))/ds
         enddo

      enddo ! 反復のループ

      ! step5 体積力の決定
      gx(:,:) = 0d0
      gy(:,:) = 0d0
      do k = 1,nb
         do j=nymin,nymax
         do i=nxmin,nxmax
            gx(i,j) = gx(i,j) + gbx(k)*w(i-xb(k))*w(j-yb(k))*dv
            gy(i,j) = gy(i,j) + gby(k)*w(i-xb(k))*w(j-yb(k))*dv
         enddo
         enddo
      enddo

      ! ===分布関数の更新 3.2
      do j=nymin,nymax
      do i=nxmin,nxmax
         do l=first+1,Last
            f(l,i,j) = f(l,i,j) + 3d0*ds*Weight(l) &
                      *(dble(ConvVelx(l))*gx(i,j) + dble(convvely(l))*gy(i,j))
         enddo
      enddo
      enddo

      

   end subroutine

   subroutine computef(f, velx, vely, dens)
      implicit none
      real(real64), intent(in)    :: f(First:Last, 1:Nx, 1:Ny)
      real(real64), intent(inout) :: velx(1:Nx, 1:Ny)
      real(real64), intent(inout) :: vely(1:Nx, 1:Ny)
      real(real64), intent(inout) :: dens(1:Nx, 1:Ny)
      integer :: i, j
      real(real64) :: f_boundary, f_exterior

      do j = 1, ny
      do i = 1, nx
         dens(i, j) = f(Center   , i, j) &
                    + f(Right    , i, j) &
                    + f(Up       , i, j) &
                    + f(Left     , i, j) &
                    + f(Down     , i, j) &
                    + f(UpRight  , i, j) &
                    + f(UpLeft   , i, j) &
                    + f(DownLeft , i, j) &
                    + f(DownRight, i, j)
      end do
      end do

      do j = 1, Ny
      do i = 1, Nx
         velx(i, j) = (f(Center   , i, j)*Convvelx(Center   ) &
                     + f(Right    , i, j)*Convvelx(Right    ) &
                     + f(Up       , i, j)*Convvelx(Up       ) &
                     + f(Left     , i, j)*Convvelx(Left     ) &
                     + f(Down     , i, j)*Convvelx(Down     ) &
                     + f(UpRight  , i, j)*Convvelx(UpRight  ) &
                     + f(UpLeft   , i, j)*Convvelx(UpLeft   ) &
                     + f(DownLeft , i, j)*Convvelx(DownLeft ) &
                     + f(DownRight, i, j)*Convvelx(DownRight))/dens(i, j)

         vely(i, j) = (f(Center   , i, j)*Convvely(Center   ) &
                     + f(Right    , i, j)*Convvely(Right    ) &
                     + f(Up       , i, j)*Convvely(Up       ) &
                     + f(Left     , i, j)*Convvely(Left     ) &
                     + f(Down     , i, j)*Convvely(Down     ) &
                     + f(UpRight  , i, j)*Convvely(UpRight  ) &
                     + f(UpLeft   , i, j)*Convvely(UpLeft   ) &
                     + f(DownLeft , i, j)*Convvely(DownLeft ) &
                     + f(DownRight, i, j)*Convvely(DownRight))/dens(i, j)

      end do
      end do

   end subroutine

   ! =====重み関数=====
   function w(r)
      real(real64)::r,w

      ! 3.5
      if(dabs(r) <= 1d0)then
         w = 0.125d0*(3d0-2d0*dabs(r)+dsqrt(1d0 + 4d0*dabs(r) - 4d0*r**2))
      else if(1d0 < dabs(r) .and. dabs(r) <= 2d0)then
         w = 0.125d0*(5d0-2d0*dabs(r)+dsqrt(-7d0 + 12d0*dabs(r) - 4d0*r**2))
      else
         w = 0d0
      endif

   end function

end program lbm
