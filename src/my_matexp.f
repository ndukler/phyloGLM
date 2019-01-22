* wrapalldmexpv:
      subroutine wrapalldmexpv(n,m,t,v,w,tol,anorm,wsp,lwsp,iwsp,
     .  liwsp,itrace,iflag,ia,ja,a,nz,res,mxstep,flag,flag2,flag3)
      implicit none
      integer,intent(inout) :: m, lwsp, liwsp, itrace, iflag
      integer,intent(inout) :: iwsp(liwsp)
      integer,intent(in) :: nz,n,ia(n+1),ja(nz), mxstep
      double precision,intent(inout) :: t, tol, anorm, w(n), v(n)
      double precision,intent(inout) :: wsp(lwsp)
      integer,intent(out) :: flag
      double precision,intent(out) :: flag2(n+1), flag3(n+1)
      double precision,intent(in) :: a(nz)
      double precision, intent(out) :: res(n*n)
      double precision ZERO, ONE
      parameter( ZERO=0.0d0, ONE=1.0d0 )
      intrinsic ABS
      integer i,j
      flag2(1) = ZERO
      flag3(1) = ZERO
      do i = 1,n
        do j = 1,n
           v(j) = ZERO
        enddo
        v(i) = ONE
        call myDMEXPV( n, m, t,v,w, tol, anorm,wsp,lwsp,
     .         iwsp,liwsp,itrace,iflag,ia,ja,a,nz, mxstep,flag)
        do j = 1,n
            res(((i-1)*n)+j) = w(j)
        enddo
        flag2(i+1)=flag
      enddo
      do i=1,n
        flag3(i+1)=flag2(i+1)+flag3(i)
      enddo 
      if(flag3(n+1).ge.1) then
        flag = 1
      else 
        flag = 0
      endif
      end subroutine wrapalldmexpv
         	  
* wrapalldgexpv:
      subroutine wrapalldgexpv(n,m,t,v,w,tol,anorm,wsp,lwsp,iwsp,
     .  liwsp,itrace,iflag,ia,ja,a,nz,res,mxstep,flag,flag2,flag3)
      implicit none

      integer,intent(inout) :: m, lwsp, liwsp, itrace, iflag
      integer,intent(inout) :: iwsp(liwsp)
      integer,intent(in) :: nz,n,ia(n+1),ja(nz), mxstep
      double precision,intent(inout) :: t, tol, anorm, v(n), w(n)
      double precision,intent(inout) :: wsp(lwsp)
      integer,intent(out) :: flag
      double precision,intent(out) :: flag2(n+1), flag3(n+1)
      double precision,intent(in) :: a(nz)
      double precision, intent(out) :: res(n*n)
      
      double precision ZERO, ONE
      parameter( ZERO=0.0d0, ONE=1.0d0 )
      intrinsic ABS
      integer i,j
      flag2(1) = ZERO
      flag3(1) = ZERO
      do i = 1,n
        do j = 1,n
          v(j) = ZERO
        enddo
        v(i) = ONE
        call myDGEXPV( n, m, t,v,w, tol, anorm,wsp,lwsp,
     .         iwsp,liwsp,itrace,iflag,ia,ja,a,nz, mxstep,flag)
        do j = 1,n
          res(((i-1)*n)+j) = w(j)
        enddo
        flag2(i+1)=flag
      enddo
      do i=1,n
        flag3(i+1)=flag2(i+1)+flag3(i)
      enddo 
      if(flag3(n+1).ge.1) then
        flag = 1
      else 
      flag = 0
      endif
      end subroutine wrapalldgexpv
        	  
* wrapdgpadm:       
      subroutine wrapdgpadm(ideg,m,t,H,ldh,wsp,lwsp,ipiv,iexph,
     .                        ns,iflag)
      implicit none

      integer,intent(inout) :: ideg,m,ldh,lwsp,iexph,ns,iflag,ipiv(m)
      double precision,intent(inout) :: t,H(ldh,m),wsp(lwsp)
      integer i,j
*      print 9001,( (H(i,j), j=1,ldh), i=1,m )
 9001 format( 4(1X,D11.2) )
 9000 format( 4(1X,I2.1) )
      call DGPADM( ideg,m,t,H,ldh,wsp,lwsp,ipiv,iexph,ns,iflag)
*      print 9001,( (wsp(iexph+(j-1)*m+i-1), j=1,m), i=1,m )
*      print 9000,( ((iexph+(j-1)*m+i-1), j=1,m), i=1,m )
      end subroutine wrapdgpadm     
        	  


