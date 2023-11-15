import multiprocessing as mp
from tqdm import tqdm, trange

def tmap(p, f, L):
    return list(tqdm(p.imap(f, L), total=len(L)))

