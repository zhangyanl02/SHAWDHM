c	Program to calculate the area ratio of each land cover type
c	Input lu.asc (1-km)
c	Output lu_class(1-9).asc (10-km)
c	by D. Yang April 2002
c
c   Variable Definitions  
c   1  = water body
c   2  = urban area
c   3  = bare land
c   4  = forest
c   5  = irrigated cropland
c   6  = nonirrigated cropland
c   7  = grassland
c   8  = shrub
c   9  = wetland
cc 10  = tundra, snow or ice  
c  
	parameter (nnc2=78,  nnr2=106)
	parameter (nnc1=780,  nnr1=1060)

	character*100, infile, outfile(9), in_dir, out_dir
	integer lu(nnr1,nnc1)
	real  ratio(10,nnr2,nnc2), tratio(nnr2,nnc2)
      DOUBLE PRECISION  x0, y0
	
      character*5 ncols,nrows
      character*12 xllcorner,yllcorner
      character*8 cellsize
      character*12 nodata
	integer nc1, nr1
c
	in_dir  =  'D:\work\dem_test/'
	out_dir =  'D:\work\parameter_test/'
c	out_dir = '/users/yang/research/hydro_model/yellow/parameter/'
c
	nc2=nnc2
	nr2=nnr2
	s2=1000  ! cellsize
c
	  infile='lu.asc'
	  call strlen(infile, ia1,ia2)
	  call strlen(in_dir, ic1,ic2)
	  open(3, file=in_dir(ic1:ic2)//infile(ia1:ia2), status='old')
          read(3, '(a,i15)') ncols,nc1
          read(3, '(a,i15)') nrows,nr1
          read(3, '(a,f16.0)') xllcorner,x0
          read(3, '(a,f16.0)') yllcorner,y0
          read(3, '(a,f16.0)') cellsize,s1
          read(3, '(a,f16.0)') nodata,znodata
          print *, 'nc1=, nr1=, x0=, y0=', nc1, nr1, x0, y0,s1,znodata	
		pause       
c
	do i=1,9
	  outfile(i) = 'lu_class'//
     $			char(48+int(i/10))//char(48+mod(i,10))//'.asc'
	  call strlen(outfile(i), ia1,ia2)
	  call strlen(out_dir, ic1,ic2)
	  open(8+i, file=out_dir(ic1:ic2)//outfile(i)(ia1:ia2),
     $             status='unknown')

	  write(8+i, '(a,i15)') ncols,nc2
        write(8+i, '(a,i15)') nrows,nr2
        write(8+i, '(a,f16.3)') xllcorner, x0
        write(8+i, '(a,f16.3)') yllcorner, y0
        write(8+i, '(a,f16.3)') cellsize,s2
        write(8+i, '(a,f16.3)') nodata,znodata
	end do


	  open(8+11, file=out_dir(ic1:ic2)//'all_class.asc',
     $              status='unknown')
	  write(8+11, '(a,i15)') ncols, nc2
        write(8+11, '(a,i15)') nrows, nr2
        write(8+11, '(a,f16.3)') xllcorner, x0
        write(8+11, '(a,f16.3)') yllcorner, y0
        write(8+11, '(a,f16.3)') cellsize, s2
        write(8+11, '(a,f16.3)') nodata, znodata


	do 200 ir2=1,nr2
	  print *,ir2
	  do ir1= 1, 10
c	  do ir1= 1, 1   ! 依实际情况而定
	    read(3,*) (lu(ir1, ic1), ic1 = 1, nnc1)
	  end do
	  do 100 ic2 = 1, nc2
	    num=0
	    ratio(1, ir2,ic2)=0.0
	    ratio(2, ir2,ic2)=0.0
	    ratio(3, ir2,ic2)=0.0
	    ratio(4, ir2,ic2)=0.0
	    ratio(5, ir2,ic2)=0.0
	    ratio(6, ir2,ic2)=0.0
	    ratio(7, ir2,ic2)=0.0
	    ratio(8, ir2,ic2)=0.0
	    ratio(9, ir2,ic2)=0.0
          ratio(10, ir2,ic2)=0.0

	    do ir1 = 1, 10  !大网格比小网格的倍数，mqh
c	    do ir1 = 1, 1
	      do ic1 = (ic2-1)*10+1, ic2*10
c	      do ic1 = (ic2-1)*1+1, ic2*1
		    if(lu(ir1,ic1) .ne. -9999) then
		      num=num+1
		      if(lu(ir1,ic1).eq.1)
     $		  ratio(1,ir2,ic2)=ratio(1,ir2,ic2)+1.0  ! type 1

		      if(lu(ir1,ic1).eq.2) 
     $		  ratio(2,ir2,ic2)=ratio(2,ir2,ic2)+1.0  ! type 2

		      if(lu(ir1,ic1).eq.3)
     $		  ratio(3,ir2,ic2)=ratio(3,ir2,ic2)+1.0  ! type 3

		      if(lu(ir1,ic1).eq.4)
     $		  ratio(4,ir2,ic2)=ratio(4,ir2,ic2)+1.0  ! type 4

		      if(lu(ir1,ic1).eq.5)
     $		  ratio(5,ir2,ic2)=ratio(5,ir2,ic2)+1.0  ! type 5

		      if(lu(ir1,ic1).eq.6)
     $		  ratio(6,ir2,ic2)=ratio(6,ir2,ic2)+1.0  ! type 6

		      if(lu(ir1,ic1).eq.7)
     $		  ratio(7,ir2,ic2)=ratio(7,ir2,ic2)+1.0  ! type 7

		      if(lu(ir1,ic1).eq.8)
     $		  ratio(8,ir2,ic2)=ratio(8,ir2,ic2)+1.0  ! type 8

		      if(lu(ir1,ic1).eq.9)
     $		  ratio(9,ir2,ic2)=ratio(9,ir2,ic2)+1.0  ! type 9
c             added by xujijun
c		      if(lu(ir1,ic1).eq.10)
c     $		  ratio(10,ir2,ic2)=ratio(10,ir2,ic2)+1.0  ! type 10

	        endif
	      end do
	    end do
c
	    do i=1,10
	      if(num .gt. 0) then
	        ratio(i,ir2,ic2) = ratio(i,ir2,ic2)/float(num)
	if(ratio(i,ir2,ic2) .gt. 1.0) print *, ratio(i,ir2,ic2), num
		    tratio(ir2,ic2)  = tratio(ir2,ic2)+ratio(i,ir2,ic2)
	      else
	        ratio(i,ir2,ic2)=znodata
		    tratio(ir2,ic2)=znodata
	      endif
	    end do
	if(tratio(ir2,ic2).ne.znodata .and. abs(tratio(ir2,ic2)-1.0).gt.
     $         0.00001) print *, tratio(ir2,ic2)
100	  continue
	  do i=1,10
	    write(8+i,*) (ratio(i,ir2,ic2), ic2=1,nc2)
	  end do
	  write(8+11,*) (tratio(ir2,ic2), ic2=1,nc2)



200	continue
	end


c
      subroutine strlen(str,l1,l2)
      character str*100
      integer i,l1,l2,k
      k=0
      do i = 1, 100
        if(k.eq.0 .and. str(i:i).NE.' ') then
          l1=i
          k=1
        elseif(k.eq.1 .and. str(i:i).EQ.' ') then
          l2 = i-1
          return
        endif
      enddo
      l2 = i
      return
      end


