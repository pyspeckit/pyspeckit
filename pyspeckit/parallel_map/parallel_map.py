"""
Parellel Map snippet by Brian Refsdal

http://www.astropython.org/snippet/2010/3/Parallel-map-using-multiprocessing
"""
import numpy
import warnings
from astropy import log
_multi=False
_ncpus=1

try:
    # May raise ImportError
    import multiprocessing
    _multi=True
   
    # May raise NotImplementedError
    _ncpus = multiprocessing.cpu_count()
except Exception as ex:
    pmap_exception = ex
    _multi=False


__all__ = ('parallel_map',)


def worker(f, ii, chunk, out_q, err_q, lock):
  """
  A worker function that maps an input function over a
  slice of the input iterable.

  :param f  : callable function that accepts argument from iterable
  :param ii  : process ID
  :param chunk: slice of input iterable
  :param out_q: thread-safe output queue
  :param err_q: thread-safe queue to populate on exception
  :param lock : thread-safe lock to protect a resource
         ( useful in extending parallel_map() )
  """
  vals = []

  # iterate over slice 
  for val in chunk:
    try:
      result = f(val)
    except Exception, e:
      err_q.put(e)
      return

    vals.append(result)

  # output the result and task ID to output queue
  out_q.put( (ii, vals) )


def run_tasks(procs, err_q, out_q, num):
  """
  A function that executes populated processes and processes
  the resultant array. Checks error queue for any exceptions.

  :param procs: list of Process objects
  :param out_q: thread-safe output queue
  :param err_q: thread-safe queue to populate on exception
  :param num : length of resultant array

  """
  # function to terminate processes that are still running.
  die = (lambda vals : [val.terminate() for val in vals
             if val.exitcode is None])

  try:
    for proc in procs:
      proc.start()

    for proc in procs:
      proc.join()

  except Exception, e:
    # kill all slave processes on ctrl-C
    die(procs)
    raise e

  if not err_q.empty():
    # kill all on any exception from any one slave
    die(procs)
    raise err_q.get()

  # Processes finish in arbitrary order. Process IDs double
  # as index in the resultant array.
  results=[None]*num;
  while not out_q.empty():
    idx, result = out_q.get()
    results[idx] = result

  try:
      # Remove extra dimension added by array_split
      return list(numpy.concatenate(results))
  except ValueError:
      return list(results)


def parallel_map(function, sequence, numcores=None):
  """
  A parallelized version of the native Python map function that
  utilizes the Python multiprocessing module to divide and 
  conquer sequence.

  parallel_map does not yet support multiple argument sequences.

  :param function: callable function that accepts argument from iterable
  :param sequence: iterable sequence 
  :param numcores: number of cores to use
  """
  if not callable(function):
    raise TypeError("input function '%s' is not callable" %
              repr(function))

  if not numpy.iterable(sequence):
    raise TypeError("input '%s' is not iterable" %
              repr(sequence))

  size = len(sequence)

  if not _multi or size == 1 or numcores == 1:
    return map(function, sequence)

  if numcores is not None and numcores > _ncpus:
      warnings.warn("Number of requested cores is greated than the "
                    "number of available CPUs.")
  elif numcores is None:
    numcores = _ncpus

  # Returns a started SyncManager object which can be used for sharing 
  # objects between processes. The returned manager object corresponds
  # to a spawned child process and has methods which will create shared
  # objects and return corresponding proxies.
  manager = multiprocessing.Manager()

  # Create FIFO queue and lock shared objects and return proxies to them.
  # The managers handles a server process that manages shared objects that
  # each slave process has access to. Bottom line -- thread-safe.
  out_q = manager.Queue()
  err_q = manager.Queue()
  lock = manager.Lock()

  # if sequence is less than numcores, only use len sequence number of 
  # processes
  if size < numcores:
    log.info("Reduced number of cores to {0}".format(size))
    numcores = size 

  # group sequence into numcores-worth of chunks
  sequence = numpy.array_split(sequence, numcores)

  procs = [multiprocessing.Process(target=worker,
           args=(function, ii, chunk, out_q, err_q, lock))
         for ii, chunk in enumerate(sequence)]

  return run_tasks(procs, err_q, out_q, numcores)


if __name__ == "__main__":
  """
  Unit test of parallel_map()

  Create an arbitrary length list of references to a single
  matrix containing random floats and compute the eigenvals
  in serial and parallel. Compare the results and timings.
  """

  import time

  numtasks = 5
  #size = (1024,1024)
  size = (512,512)

  vals = numpy.random.rand(*size)
  f = numpy.linalg.eigvals

  iterable = [vals]*numtasks

  print ('Running numpy.linalg.eigvals %iX on matrix size [%i,%i]' %
      (numtasks,size[0],size[1]))

  tt = time.time()
  presult = parallel_map(f, iterable)
  print 'parallel map in %g secs' % (time.time()-tt)

  tt = time.time()
  result = map(f, iterable)
  print 'serial map in %g secs' % (time.time()-tt)

  assert (numpy.asarray(result) == numpy.asarray(presult)).all()
