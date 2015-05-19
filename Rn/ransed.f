      subroutine ransed(iarg,lu,iseed,iflag)
c
c  routine to set and store seed ix for ranv.
c  iarg=-2  write ix onto tape1
c  iarg=-1  read ix from tape 1
c  iarg=0  use default value (11261949)
c  iarg.gt.0  ix=iarg
c       values are written on logical unit lu
c
      common/rancom/ix
c
      if(iarg.lt.-1) go to 200
c  setting the seed
      if(iarg.eq.0) go to 100
      ix=iarg
      if(iarg.gt.0) go to 100
      if(iflag.eq.0)then
        read(1,1)ix
      else
	ix = iseed
      endif
100   continue
      !write(lu,97)ix
      return
c
c  writing ix on tape1
200   rewind (1)
      write(1,1)ix
      !write(lu,98)ix
      return
c
1      format(1x,i12)
97      format(1x,25hinitial value of seed is ,i12)
98      format(/1x,23hfinal value of seed is ,i12)
c
      end
