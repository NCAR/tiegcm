#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>
#include <stdlib.h>
#include <errno.h>
#include <signal.h>

#include <sys/types.h>
#include <sys/wait.h>
#include <sys/poll.h>

#define MAX_LABEL_LEN 256

char label[MAX_LABEL_LEN];


void do_io(int rfd, int wfd, int*newline, int *eof)
{
  char rbuf[1024];
  char obuf[1024+MAX_LABEL_LEN];
  int n;
  int err;
  do {
    n = read(rfd,rbuf,sizeof(rbuf));
    err = errno;
    if (n > 0) {
      /* copy read buffer,  inserts label after every new line */
      int i,o;
      for (i =0,o=0;i<n;i++,o++) {
	if (o >= sizeof(obuf) - MAX_LABEL_LEN - 1) {
	  write(wfd,obuf,o);
	  o = 0;
	}
	if (*newline) {
	  o += sprintf(obuf+o,"%s",label);
	  *newline = 0;
	}
	obuf[o] = rbuf[i];
	if (rbuf[i] == '\n') {
	  *newline = 1;
	}
      }
      write(wfd,obuf,o);
    }
  } while (n > 0);
  *eof = (n == 0) || (n == -1 && err != EAGAIN);
}

int main(int argc,char *argv[])
{

  pid_t child;
  int status;
  int pipe_stdout[2];
  int pipe_stderr[2];
  int cc;
  struct pollfd ufds[4];
  int stdout_nl = 1;
  int stderr_nl = 1;
  int stdout_eof = 0;
  int stderr_eof = 0;
  char *s;
  if ((s = getenv("GMPI_ID"))) {
    sprintf(label,"%s:",s);
  }
  signal(SIGPIPE, SIG_IGN);
  /* signal(SIGCHLD, handler); */
  if (pipe(pipe_stdout) != 0 || pipe(pipe_stderr) != 0) {
    perror("pipe");
    exit(1);
  }
  child = fork();
  if (child == 0) {
    close(pipe_stdout[0]);
    close(pipe_stderr[0]);
    dup2(pipe_stdout[1],1);
    dup2(pipe_stderr[1],2);
    execvp(argv[1],argv+1);
    _exit(2);
  } else if (child < 0) {
    perror("fork");
    _exit(3);
  }
  close(pipe_stdout[1]);
  close(pipe_stderr[1]);
  fcntl(pipe_stdout[0], F_SETFL, fcntl(pipe_stdout[0],F_GETFL) | O_NONBLOCK);
  fcntl(pipe_stderr[0], F_SETFL, fcntl(pipe_stderr[0],F_GETFL) | O_NONBLOCK);
  ufds[2].events = 0;
  do {
    int n = 0;
    int stdout_idx, stderr_idx;
    if (!stdout_eof) {
      ufds[n].fd = pipe_stdout[0];
      ufds[n].events = POLLIN;
      stdout_idx = n;
      n += 1;
    }
    if (!stderr_eof) {
      ufds[n].fd = pipe_stderr[0];
      ufds[n].events = POLLIN;
      stderr_idx = n;
      n += 1;
    }
    ufds[n].fd = 1;
    ufds[n].events = 0;
    n += 1;
    ufds[n].fd = 2;
    ufds[n].events = 0;
    n += 1;
    cc = poll(ufds,n,-1);
    if (cc < 0 && errno == EINTR) {
      continue;
    } else if (cc < 0) {
      perror("poll");
      kill(child,SIGTERM);
      exit(1);
    }
    if (!stdout_eof && ufds[stdout_idx].revents) {
      do_io(pipe_stdout[0],1,&stdout_nl, &stdout_eof);
    }
    if (!stderr_eof && ufds[stderr_idx].revents) {
      do_io(pipe_stderr[0],2,&stderr_nl, &stderr_eof);
    }
    if (ufds[n - 2].revents) {
      break;
    }
    if (ufds[n - 1].revents) {
      break;
    }
  } while (stdout_eof + stderr_eof < 2);
  kill(child,SIGTERM);
  waitpid(child,&status,0);
  if (WIFEXITED(status)) {
    if (WEXITSTATUS(status)) {
      fprintf(stderr,"%s%s:exit(%d)\n",label, argv[1], WEXITSTATUS(status));
    }
    return WEXITSTATUS(status);
  } else {
    fprintf(stderr,"%s%s:signal %d\n", label, argv[1], WTERMSIG(status));
    raise(WTERMSIG(status));
  }
  _exit(4);
}
