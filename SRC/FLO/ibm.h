#if PARALLEL

#define PARALLEL_START() CALL MPI_INIT(err_code)
#define PARALLEL_FINISH()CALL MPI_FINALIZE(err_code)  
#define GET_NPROCS(comm,psize) CALL MPI_COMM_SIZE(comm,psize,err_code)
#define GET_PROCESSOR_RANK(comm,rank) CALL MPI_COMM_RANK(comm, rank, err_code)
#define BROADCAST_INT(item,nitems,from,comm)   CALL MPI_BCAST(item,nitems,MPI_INTEGER,from,comm,err_code)
#define BROADCAST_DOUBLE(item,nitems,from,comm)   CALL MPI_BCAST(item,nitems,MPI_DOUBLE_PRECISION,from,comm,err_code)
#define BROADCAST_REAL(item,nitems,from,comm)   CALL MPI_BCAST(item,nitems,MPI_DOUBLE_REAL,from,comm,err_code)
#define BROADCAST_CHARARR(str,strlen,from,comm)   CALL MPI_BCAST(str,strlen,MPI_CHARACTER,from,comm,err_code)
#define BROADCAST_LOGICAL(item,nitems,from,comm)   CALL MPI_BCAST(item,nitems,MPI_LOGICAL,from,comm,err_code)
#define BROADCAST_STRING(str,strlen,from,comm)				\
  if(I_AM_NODE_ZERO) strlen = LEN_TRIM(str); BROADCAST_INT(strlen,1,from,comm); BROADCAST_CHARARR(str,strlen,from,comm)
#define SEND_STRING(str,strlen,to,sitag,sctag,comm) strlen = LEN_TRIM(str); CALL MPI_SEND(strlen,1,MPI_INTEGER,to,sitag,comm,err_code); CALL MPI_SEND(str,strlen,MPI_CHARACTER,to,sctag,comm,err_code)
#define RECV_STRING(str,strlen,from,ritag,rctag,comm,status) CALL MPI_RECV(strlen,1,MPI_INTEGER,from,ritag,comm,status,err_code); CALL MPI_RECV(str,strlen,MPI_CHARACTER,from,rctag,comm,status,err_code)

#define SEND_CHARACTER(sitem,nsitems,to,stag,comm) CALL MPI_SEND(sitem,nsitems,MPI_CHARACTER, to,stag,comm,err_code)
#define RECV_CHARACTER(ritem,nritems,from,rtag,comm,status) CALL MPI_RECV(ritem,nritems, MPI_CHARACTER,from,rtag,comm,status,err_code)

#define SEND_INT(sitem,nsitems,to,stag,comm) CALL MPI_SEND(sitem,nsitems,MPI_INTEGER, to,stag,comm,err_code)
#define RECV_INT(ritem,nritems,from,rtag,comm,status) CALL MPI_RECV(ritem,nritems, MPI_INTEGER,from,rtag,comm,status,err_code)

#define CREATE_CART_TOPOLOGY(commold,ndims,psize,isperiodic,reorder,newcomm) CALL MPI_CART_CREATE(commold,ndims,psize,isperiodic,reorder,newcomm,err_code)
#define GET_SHIFT_PROCS(comm,dir,shift,src,dest) CALL MPI_CART_SHIFT(comm,dir,shift,src,dest,err_code)
#define GLOBAL_INT_SUM(sendbuf,recvbuf,nitems,comm) CALL MPI_ALLREDUCE(sendbuf,recvbuf,nitems,MPI_INTEGER,MPI_SUM,comm,err_code)
#define GLOBAL_LONGINT_SUM(sendbuf,recvbuf,nitems,comm) CALL MPI_ALLREDUCE(sendbuf,recvbuf,nitems,MPI_INTEGER8,MPI_SUM,comm,err_code)
#define GLOBAL_DOUBLE_SUM(sendbuf,recvbuf,nitems,comm) CALL MPI_ALLREDUCE(sendbuf,recvbuf,nitems,MPI_DOUBLE_PRECISION,MPI_SUM,comm,err_code)
#define GLOBAL_COMPLEX_SUM(sendbuf,recvbuf,nitems,comm) CALL MPI_ALLREDUCE(sendbuf,recvbuf,nitems,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,err_code)
#define GLOBAL_DOUBLE_MAX(sendbuf,recvbuf,nitems,comm) CALL MPI_ALLREDUCE(sendbuf,recvbuf,nitems,MPI_DOUBLE_PRECISION,MPI_MAX,comm,err_code)
#define GLOBAL_DOUBLE_MIN(sendbuf,recvbuf,nitems,comm) CALL MPI_ALLREDUCE(sendbuf,recvbuf,nitems,MPI_DOUBLE_PRECISION,MPI_MIN,comm,err_code)

#define RSENDRECV(sitem,nsitems,to,stag,ritem,nritems,from,rtag,comm,status) CALL MPI_SENDRECV(sitem,nsitems,MPI_DOUBLE_PRECISION,to,stag,ritem,nritems,MPI_DOUBLE_PRECISION,from,rtag,comm,status,err_code)
#define CSENDRECV(sitem,nsitems,to,stag,ritem,nritems,from,rtag,comm,status) CALL MPI_SENDRECV(sitem,nsitems,MPI_DOUBLE_COMPLEX,to,stag,ritem,nritems,MPI_DOUBLE_COMPLEX,from,rtag,comm,status,err_code)

#define CREATE_XY_LSLICE(count,newtype) CALL MPI_TYPE_CONTIGUOUS(count,MPI_LOGICAL,newtype,err_code)
#define CREATE_XY_RSLICE(count,newtype) CALL MPI_TYPE_CONTIGUOUS(count,MPI_DOUBLE_PRECISION,newtype,err_code)

#define CREATE_XZ_LSLICE(count,blocklength,stride,newtype) CALL MPI_TYPE_VECTOR(count,blocklength,stride,MPI_LOGICAL,newtype,err_code)
#define CREATE_XZ_RSLICE(count,blocklength,stride,newtype) CALL MPI_TYPE_VECTOR(count,blocklength,stride,MPI_DOUBLE_PRECISION,newtype,err_code)

#define COMMIT(dtype) CALL MPI_TYPE_COMMIT(dtype,err_code)
#define VECSENDRECV(sitem,nsitems,vectype,to,stag,ritem,nritems,from,rtag,comm,status) CALL MPI_SENDRECV(sitem,nsitems,vectype,to,stag,ritem,nritems,vectype,from,rtag,comm,status,err_code)
#define VECSENDRECV2(sitem,nsitems,vectype1,to,stag,ritem,nritems,vectype2,from,rtag,comm,status) CALL MPI_SENDRECV(sitem,nsitems,vectype1,to,stag,ritem,nritems,vectype2,from,rtag,comm,status,err_code)
#define BARRIER(comm) CALL MPI_BARRIER(comm,err_code)

#else

#define PARALLEL_START() 
#define PARALLEL_FINISH() 
#define GET_NPROCS(comm,psize) nproc = 1
#define GET_PROCESSOR_RANK(comm,rank) rank = 0
#define BROADCAST_INT(item,nitems,from,comm)
#define BROADCAST_DOUBLE(item,nitems,from,comm)
#define BROADCAST_REAL(item,nitems,from,comm)
#define BROADCAST_CHARARR(str,strlen,from,comm)
#define BROADCAST_STRING(str,strlen,from,comm)
#define SEND_STRING(str,strlen,to,sitag,sctag,comm)		
#define RECV_STRING(str,strlen,from,ritag,rctag,comm,status)		
#define BROADCAST_LOGICAL(item,nitems,from,comm)
#define CREATE_CART_TOPOLOGY(commold,ndims,psize,isperiodic,reorder,newcomm)
#define GET_SHIFT_PROCS(comm,dir,shift,src,dest)
#define GLOBAL_INT_SUM(sendbuf,recvbuf,nitems,comm) recvbuf = sendbuf
#define GLOBAL_LONGINT_SUM(sendbuf,recvbuf,nitems,comm) recvbuf = sendbuf
#define GLOBAL_DOUBLE_MAX(sendbuf,recvbuf,nitems,comm) recvbuf = sendbuf
#define GLOBAL_DOUBLE_MIN(sendbuf,recvbuf,nitems,comm) recvbuf = sendbuf
#define GLOBAL_DOUBLE_SUM(sendbuf,recvbuf,nitems,comm) recvbuf = sendbuf
#define GLOBAL_COMPLEX_SUM(sendbuf,recvbuf,nitems,comm) recvbuf = sendbuf
#define RSENDRECV(sitem,nsitems,to,stag,ritem,nritems,from,rtag,comm,status)
#define CSENDRECV(sitem,nsitems,to,stag,ritem,nritems,from,rtag,comm,status)
#define VECSENDRECV(sitem,nsitems,vectype,to,stag,ritem,nritems,from,rtag,comm,status)
#define VECSENDRECV2(sitem,nsitems,vectype1,to,stag,ritem,nritems,vectype2,from,rtag,comm,status)
#define CREATE_XY_LSLICE(count,newtype)
#define CREATE_XY_RSLICE(count,newtype)
#define CREATE_XZ_RSLICE(count,blocklength,stride,newtype)
#define CREATE_XZ_LSLICE(count,blocklength,stride,newtype)
#define SEND_CHARACTER(sitem,nsitems,to,stag,comm)
#define RECV_CHARACTER(ritem,nritems,from,rtag,comm,status)
#define SEND_INT(sitem,nsitems,to,stag,comm)
#define RECV_INT(ritem,nritems,from,rtag,comm,status)
#define BARRIER(comm)
#define COMMIT(newtype)
#endif
#define I_AM_NODE_ZERO (myid.eq.0)


