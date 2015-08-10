#!/home/dence/applications/anaconda/bin/python

from multiprocessing import Pool

def foo(args):
	print args[0] + args[1]

if __name__ == '__main__':

	pool = Pool(10)
	my_ints=[(1,2),(3,4),(5,6)]
	pool.map(foo, my_ints)
