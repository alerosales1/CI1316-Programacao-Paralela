/********************************************************************************/
/*                                                                              */
/* A leitura dos querys é feita sequêncialmente, mas a busca em cada sequência  */
/* é efetuada em paralela e distribuida entre os threads.                       */
/*                                                                              */
/* versão otimizada                                                             */
/*                                                                              */
/********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <math.h>
#include <float.h>
#include <sys/timeb.h>

#define NUM_THREADS 4 //defino manualmente, #pragma opm parallel (NUM_THREADS)
#define MAX 256 

typedef struct Sequencia
 {
  char *desc; //ponteiro descrição
  char *seq;  // ponteiro sequencia 
  int  lenght; //tamanho 
 } sequencia;

struct timeb start,stop;

FILE *fdatabase, *fquery, *fout; 
char *bases;  
int tab[1000];
int num;
sequencia *seq = NULL;  
sequencia *tab_seq[1000]; // tabela de sequencia que guarda a informação dna.in

// Boyers-Moore-Hospool-Sunday algorithm for string matching
int bmhs(char *string, int n, char *substr, int m)
 {
  int d[MAX];
  int i, j, k;

  // pre-processing
  for ( j = 0; j < MAX; j++) d[j] = m + 1;
  
  for (j = 0; j < m; j++) d[(int) substr[j]] = m - j;

  // searching
  i = m - 1;

  while ( i < n )
   {
    k = i;
    j = m - 1;

    while ( (j >= 0) && (string[k] == substr[j]) ) { j--; k--; }

    if ( j < 0 ) return k + 1;

    i = i + d[(int) string[i + 1]];
   }

  return -1;
 }


void openfiles()
 {
  fdatabase = fopen("dna.in","r+");
  //fdatabase = fopen("dmel-all-transcript-r5.25.fasta","r+");

  if ( fdatabase == NULL )
   {
    perror("dna.in");
    exit(EXIT_FAILURE);
   }

  fquery = fopen("query.in","r");

  if ( fquery == NULL )
   {
    perror("query.in");
    exit(EXIT_FAILURE);
   }

  fout = fopen("dna.out","w");

  if (fout == NULL)
   {
    perror("dna.out");
    exit(EXIT_FAILURE);
   }
 }

void closefiles()
 {
  fflush(fdatabase);
  fclose(fdatabase);

  fflush(fquery);
  fclose(fquery);

  fflush(fout);
  fclose(fout);
 }

static inline void remove_eol(char *line)
 {
  int i = strlen(line) - 1;

  while ( line[i] == '\n' || line[i] == '\r' )
   {
    line[i] = 0;
    i--;
   }
 }
 
int index_base() // leitura e indexação das sequencias de DNA
 {
  char line[100];
  int i,j,ind;
  int num;

  ind = 0;
  j = 1;
  num = 1;

  tab[0] = 0;
  fseek(fdatabase,0,SEEK_SET);
  fgets(line,100,fdatabase);

  ind += strlen(line);

  while ( !feof(fdatabase) )
   {
    fgets(line,100,fdatabase);
    ind += strlen(line);

    do
     {
      if ( fgets(line, 100, fdatabase) == NULL ) break;
      
      if ( line[0] != '>' ) ind += strlen(line); //inicio de cada sequencia 
      else num++;
     }
    while (line[0] != '>'); //final da sequencia 

    tab[j] = ind;
    j++;
    ind += strlen(line);
   }
    
  tab[j] = ind;
 
  return num; // número de sequência de DNA no arquivo
 }

void prepare_base() // armazamento e preparação das sequencias para a função bmhs
 {
  int ind,lenght,quot,j,k;
  char *ptr,*seq_dna;
  char buf[25000];

  bases = (char*) malloc(sizeof(char) * 1000001);

  if ( bases == NULL )
   {
    printf("erro de alocação de memória : bases\n");
    exit (-1);
   }

  for ( ind = 0 ; ind < num ; ind++ ) // para todas as sequencias
   { 
    lenght = tab[ind+1]-tab[ind]; // comprimento da sequencia
    ptr = bases + ((int) floor ( 1000001.0 / num ) * ind); // partilha da memoria entre as sequencias

    fseek(fdatabase,tab[ind],SEEK_SET);       // posicionamento no inicio da sequencia SEEK_SET seta posição no TAB
    fread(ptr,sizeof(char),lenght,fdatabase); // leitura da sequencia

    // marca o final da sequência
    *(ptr+lenght) = '\0';
         
    // eliminação do \r\n no final da sequência
    lenght -= 2;
    *(ptr+lenght) = '\0';

    // identificação do final da descrição da sequência
    for ( j = 81 ; *(ptr+j) != '\n' ; j--);
    *(ptr+j-1) = '\0';

    tab_seq[ind] = (sequencia *) malloc(sizeof(sequencia)); // alocação memoria para o elemento
    tab_seq[ind]->desc = ptr;     // descrição da sequencia       
    tab_seq[ind]->seq  = ptr+j+1; // inicio da sequência
    seq_dna = tab_seq[ind]->seq;

    // determinação da posição do \r\n na sequência
    lenght = strlen(seq_dna);           // comprimento da sequencia
    quot = (int) floor (lenght/82.0);   // grupos de oitenta carateres + \r\n

    // eliminação dos \r\n na sequência
    for ( k = quot-1 ; k >= 0 ; k-- ) 
     {
      strcpy(buf,seq_dna+82+82*k);
      *(seq_dna+80+82*k) = '\0';
      strcat(seq_dna,buf); 
     }

    tab_seq[ind]->lenght = strlen(seq_dna);  // guarda o comprimento da sequencia
   } // end for ( ind = 0 ; ind < num ; ind++ )
 }

int main (void)
 {
  char *desc_dna;
  char *seq_dna; 
  char *seq_query;
  char desc_query[100];
  char line[100];
  int  i,j,k,found,result,ind,id,n,query_lenght;
  char ch,resultado[25000],buf[25000];

  ftime(&start); // marcar o tempo inicial

  seq_query = (char*) malloc(sizeof(char) * 1000001); // alocação memoria para a sequencia de busca

  if ( seq_query == NULL )
   {
    perror("malloc seq_query");
    exit(EXIT_FAILURE);
   }

  openfiles();

  num = index_base(); // indexação da base de sequencias
  prepare_base();     // preparação da base para a função bmhs()

  n = (int) floor(((double) num) / ((double) NUM_THREADS));

  fgets(desc_query,100,fquery); // leitura da descrição da sequencia de busca
  remove_eol(desc_query);

  while ( !feof(fquery) )
   {
    fprintf(fout,"\n%s\n",desc_query);
    resultado[0] = '\0';

    // read query string
    fgets(line,100,fquery);
    remove_eol(line);

    seq_query[0] = 0;
    i = 0;

    do
     {
      strcat(seq_query+i,line);

      if ( fgets(line,100,fquery) == NULL ) break;
      remove_eol(line);

      i += 80;
     }
    while ( line[0] != '>' );

    strcpy(desc_query,line);
    query_lenght = strlen(seq_query);
    found = 0;

    #pragma omp parallel num_threads(NUM_THREADS) private(ind,id,i,result) shared(n,num,seq_query,query_lenght,found,tab_seq)
     {
      id = omp_get_thread_num();  // identificator do thread

      for ( i = 0 ; i <= n ; i++ ) // A busca nas sequências são partilhadas entre os threads
       { 
        ind = i*NUM_THREADS + id; // escolha da sequencia a ser analisada pelo threat, qual ser tratada

        if ( ind < num )
         {
          result = bmhs(tab_seq[ind]->seq,tab_seq[ind]->lenght,seq_query,query_lenght); //resultado da sequencia bmhs

          if (result > 0)
           {
            sprintf(buf,"%s\n%d\n",tab_seq[ind]->desc,result); // salvar o resultado da busca
            strcat(resultado,buf);
            found++;
           } 
         } // end if ( ind < num )
       } // end for ( i = 0 ; i <= n ; i++ ) 
     } // end #pragma omp parallel

    // impressão do resultado da busca
    if (found) fprintf(fout,"%s",resultado);
    else       fprintf(fout,"NOT FOUND\n");
   } // end while ( !feof(fquery) )
    
  closefiles();
  free(seq_query);
  free(bases);
  for ( i = 0 ; i < ind ; i++ ) free(tab_seq[i]);

  ftime(&stop); // marcar o tempo final
  printf("tempo de execução: %.3f s\n",stop.time-start.time + (stop.millitm-start.millitm)/1000.0);

  return EXIT_SUCCESS;
 }

