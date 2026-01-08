#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int
score(const char * seq, int seqlen, const int * posx, const int * posy)
{
  int pairs = 0;
  for (int ii = 0; ii < seqlen; ++ ii) {
    if (seq[ii] == 'H') {
      const int ix = posx[ii];
      const int iy = posy[ii];
      for (int jj = 0; jj < seqlen; ++ jj) {
	if ((jj < ii - 1 || ii + 1 < jj) && seq[jj] == 'H') {
	  const int jx = posx[jj];
	  const int jy = posy[jj];
	  if ((jx == ix && jy == iy + 1) || (jy == iy && jx == ix + 1)) {
	    ++ pairs;
	  }
	}
      }
    }
  }
  return pairs;
}

int
trydir(int len, char * dir, int * posx, int * posy, char try, int dx, int dy)
{
  const int newx = posx[len - 1] + dx;
  const int newy = posy[len - 1] + dy;

  for (int ii = len - 1; ii >= 0; -- ii) {
    if (posx[ii] == newx && posy[ii] == newy) {
      return 0;
    }
  }

  dir[len] = try;
  posx[len] = newx;
  posy[len] = newy;

  return 1;
}

void
enumerate(const char * seq, int seqlen, int len, char * dir, int * posx, int * posy, int is_first, int is_flat, unsigned int * count)
{
  if (len == seqlen) {
    ++ (*count);
    const int sc = score(seq, seqlen, posx, posy);
    if (sc > 0) {
      printf("%2d ", sc);
      for (int ii = 1; ii < len; ++ ii) {
	      printf("%c", dir[ii]);
      }
      printf("\n");
    }
  } else {
    if (trydir(len, dir, posx, posy, 'R', 1, 0)) {
      enumerate(seq, seqlen, len + 1, dir, posx, posy, 0, is_flat, count);
    }
    if (!is_first) {
      if (trydir(len, dir, posx, posy, 'U', 0, 1)) {
	enumerate(seq, seqlen, len + 1, dir, posx, posy, 0, 0, count);
      }
      if (!is_flat) {
	if (trydir(len, dir, posx, posy, 'L', -1, 0)) {
	  enumerate(seq, seqlen, len + 1, dir, posx, posy, 0, 0, count);
	}
	if (trydir(len, dir, posx, posy, 'D', 0, -1)) {
	  enumerate(seq, seqlen, len + 1, dir, posx, posy, 0, 0, count);
	}
      }
    }
  }
}

int
main(int argc, const char ** argv)
{
  if (argc != 2) {
    return -1;
  }

  const char * seq = argv[1];
  const int seqlen = strlen(seq);

  for (int ii = 0; ii < seqlen; ++ ii) {
    if (seq[ii] != 'H' && seq[ii] != 'P') {
      return -1;
    }
  }

  char * dir = (char *)malloc(seqlen * sizeof(char));
  int * posx = (int *)malloc(seqlen * sizeof(int));
  int * posy = (int *)malloc(seqlen * sizeof(int));

  dir[0] = 'X';
  posx[0] = 0;
  posy[0] = 0;

  unsigned int count = 0;

  enumerate(seq, seqlen, 1, dir, posx, posy, 1, 1, &count);

  // printf("%u\n", count);

  free(posy);
  free(posx);
  free(dir);
}
