#include <R.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>

#define MAX_BARCODE 1000
#define MAX_HAIRPIN 10000
#define SEQ_LEN 100
#define BLOCKSIZE 10000000

typedef struct {
   char   *sequence;
   int    original_pos;
} a_barcode;

typedef struct {
   char   *sequence;
   int    original_pos;
   long   count;
} a_hairpin;

a_barcode *barcodes[MAX_BARCODE];
a_hairpin *hairpins[MAX_HAIRPIN];

int num_barcode;
int num_hairpin;
long num_read;
long summary[MAX_HAIRPIN][MAX_BARCODE];
int barcode_start;
int barcode_end;
int barcode_length;
int hairpin_start;
int hairpin_end;
int hairpin_length;
int allow_shifting;
int shifting_n_base; 
int allow_mismatch;
int num_mismatch_hairpin;
int barcode_n_mismatch;
int hairpin_n_mismatch;
int allow_shifted_mismatch;
int isverbose;

long barcodecount;
long hairpincount;
long bchpcount;

a_hairpin *mismatch_hairpins[MAX_HAIRPIN];
int *barcodeindex;
int *hairpinindex;

void
Read_In_Barcodes(char* filename){
  FILE *fin;
  char * line = NULL;
  size_t len = 1000;
  char *readline;

  fin = fopen(filename,"r");
  line = (char *)malloc(len+1);
  a_barcode *new_barcode;

  int count = 0;
  while ((readline = fgets(line, len, fin)) != NULL){
    count++;
    new_barcode = (a_barcode *)malloc(sizeof(a_barcode));
    new_barcode->sequence = (char *)malloc(SEQ_LEN * sizeof(char));

    new_barcode->original_pos = count;
    strncpy(new_barcode->sequence, line, barcode_length);
    barcodes[count] = new_barcode;
  }
  fclose(fin);
  num_barcode = count;
  free(line);
  Rprintf(" -- Number of Barcodes is : %d\n", num_barcode);
}

int
locate_barcode(char *a_barcode){
  int imin, imax, imid;
  imin = 1;
  imax = num_barcode;

  while (imax >= imin) {
    imid = (imax + imin) / 2;
    
    if (strncmp(barcodes[imid]->sequence, a_barcode, barcode_length) < 0)
      imin = imid + 1;
    else if (strncmp(barcodes[imid]->sequence, a_barcode, barcode_length) > 0)
      imax = imid - 1;
    else
      return barcodes[imid]->original_pos;
  }
  return -1;  
}

int
locate_hairpin(char *a_hairpin){
  int imin, imax, imid;
  imin = 1;
  imax = num_hairpin;

  while (imax >= imin) {
    imid = (imax + imin) / 2;
    
    if (strncmp(hairpins[imid]->sequence, a_hairpin, hairpin_length) < 0)
      imin = imid + 1;
    else if (strncmp(hairpins[imid]->sequence, a_hairpin, hairpin_length) > 0)
      imax = imid - 1;
    else
      return hairpins[imid]->original_pos;
  }
  return -1;  
}

int
Valid_Match(char *sequence1, char *sequence2, int length, int threshold){
  int i_base;
  int mismatchbasecount = 0;
  for (i_base = 0; i_base < length; i_base++) {
    if (sequence1[i_base] != sequence2[i_base])
      mismatchbasecount++;		  
  }
  if (mismatchbasecount <= threshold)
    return 1;
  else
    return -1;
}

int
locate_mismatch_barcode(char *a_barcode){
  int i;
  int match_index = -1;
  for (i = 1; i <= num_barcode; i++){
    if (Valid_Match(a_barcode, barcodes[i]->sequence, barcode_length, barcode_n_mismatch) > 0) {
      match_index = barcodes[i]->original_pos;
      break;
    }
  }
  return match_index;
}

int
locate_mismatch_hairpin(char *a_hairpin){
  int i;
  int match_index = -1;
  for (i = 1; i <= num_mismatch_hairpin; i++){
    if (Valid_Match(a_hairpin, mismatch_hairpins[i]->sequence, hairpin_length, hairpin_n_mismatch) > 0) {
      match_index = mismatch_hairpins[i]->original_pos;
      break;
    }
  }
  return match_index;
}


void
Sort_Barcodes(void){
  int i, j;
  a_barcode *temp;
  for(i = 1; i < num_barcode; i++){
    for(j = i+1; j <= num_barcode; j++){
      if (strcmp(barcodes[i]->sequence, barcodes[j]->sequence) > 0){
	temp = barcodes[i];
	barcodes[i] = barcodes[j];
	barcodes[j] = temp;
      }
    }
  }
}

void
Read_In_Hairpins(char *filename){
  FILE *fin;
  char * line = NULL;
  size_t len = 1000;
  char *readline;

  fin = fopen(filename,"r");
  line = (char *)malloc(len+1);
  a_hairpin *new_hairpin;

  int count = 0;
  while ((readline = fgets(line, len, fin)) != NULL){
    count++;
    new_hairpin = (a_hairpin *)malloc(sizeof(a_hairpin));
    new_hairpin->sequence = (char *)malloc(SEQ_LEN * sizeof(char));
    new_hairpin->original_pos = count;
    new_hairpin->count = 0;
    strncpy(new_hairpin->sequence, line, hairpin_length);
    hairpins[count] = new_hairpin;
  }
  fclose(fin);
  num_hairpin = count;
  free(line);
  Rprintf(" -- Number of Hairpin is : %d\n", num_hairpin);
}

void
Sort_Hairpins(void){
  int i, j;
  a_hairpin *temp;
  for(i = 1; i < num_hairpin; i++){
    for(j = i+1; j <= num_hairpin; j++){
      if (strcmp(hairpins[i]->sequence, hairpins[j]->sequence) > 0){
	temp = hairpins[i];
	hairpins[i] = hairpins[j];
	hairpins[j] = temp;
      }
    }
  }
}

long Count_Reads(char *filename) {
  FILE *freads = fopen(filename, "r");
  char * line = NULL;
  char *readline;
  size_t len = 1000;
  line = (char *)malloc(sizeof(char) * (len + 1));
  if(freads == NULL){
    fclose(freads);
    return 0;
  }
  long line_count = 0;
  while ((readline = fgets(line, len, freads)) != NULL){
    line_count++;
  }
  fclose(freads);
  free(line);
  return line_count / 4;
}


void
Process_Hairpin_Reads(char *filename){
  FILE *fin;
  char * line = NULL;
  size_t len = 1000;
  char *readline;
  long num_read_thisfile = 0;

  fin = fopen(filename,"r");
  line = (char *)malloc(sizeof(char) * (len+1));
  
  if (isverbose)
    Rprintf("Processing reads in %s.\n", filename);

  char * this_barcode;
  char * this_hairpin;
  this_barcode = (char *)malloc(SEQ_LEN * sizeof(char));
  this_hairpin = (char *)malloc(SEQ_LEN * sizeof(char));

  long line_count = 0;

  int barcode_index;
  int hairpin_index;

  while ((readline = fgets(line, len, fin)) != NULL){
    line_count++;  
    if ((line_count % 4) != 2)
      continue;
    
    if ((isverbose) && (num_read_thisfile % BLOCKSIZE == 0))
      Rprintf(" -- Processing %d million reads\n", (num_read_thisfile / BLOCKSIZE + 1) * 10);
    num_read++;
    num_read_thisfile++;
    strncpy(this_barcode, line + barcode_start - 1, barcode_length);
    barcode_index = locate_barcode(this_barcode);

    strncpy(this_hairpin, line + hairpin_start - 1, hairpin_length);
    hairpin_index = locate_hairpin(this_hairpin);

    if ((hairpin_index <= 0) && (allow_shifting > 0)){
      // Check if given hairpin can be mapped to a shifted location. 
      int index;
      // check shifting leftwards
      for (index = 1; index <= shifting_n_base; index++){
	strncpy(this_hairpin, line + hairpin_start - 1 - index, hairpin_length);
	hairpin_index = locate_hairpin(this_hairpin);
	if (hairpin_index > 0)
	  break;
      }
      // check shifting rightwards
      if (hairpin_index <= 0){
	for (index = 1; index <= shifting_n_base; index++){
	  strncpy(this_hairpin, line + hairpin_start - 1 + index, hairpin_length);
	  hairpin_index = locate_hairpin(this_hairpin);
	  if (hairpin_index > 0)
	    break;
	}
      }
    }

    if (barcode_index > 0)
      barcodecount++;

    if (hairpin_index > 0){
      hairpincount++;
      hairpins[hairpin_index]->count++;
    }
       
    if ((barcode_index > 0) && (hairpin_index > 0)) {
      summary[hairpin_index][barcode_index]++;
      bchpcount++;
    }
    
    barcodeindex[num_read] = barcode_index;
    hairpinindex[num_read] = hairpin_index;
  }
  if (isverbose)
    Rprintf("Number of Read in file %s is : %ld\n", filename, num_read_thisfile);  
  fclose(fin);
  free(line);
  free(this_barcode);
  free(this_hairpin);  
}

void
Create_Mismatch_Hairpins_List(void){
  int i;
  num_mismatch_hairpin = 0;

  for (i = 1; i <= num_hairpin; i++){
    if (hairpins[i]->count == 0){
      num_mismatch_hairpin++;
      mismatch_hairpins[num_mismatch_hairpin] = hairpins[i];
    }
  }
  Rprintf("\nThere are %d hairpins without exact sequence match.\n", num_mismatch_hairpin);
}


void
Process_Mismatch(char *filename){

  FILE *fin;
  char * line = NULL;
  size_t len = 1000;
  char *readline;
  long num_read_thisfile = 0;

  fin = fopen(filename,"r");
  line = (char *)malloc(len+1);

  if (isverbose)
    Rprintf("Processing reads in %s, considering sequence mismatch. \n", filename);

  char * this_hairpin;
  char * this_barcode;
  this_hairpin = (char *)malloc(SEQ_LEN * sizeof(char));
  this_barcode = (char *)malloc(SEQ_LEN * sizeof(char));

  long line_count = 0;

  int new_barcode_index;
  int new_hairpin_index;
  while ((readline = fgets(line, len, fin)) != NULL){
    line_count++;
    if ((line_count % 4) != 2)
      continue;

    if ((isverbose) && (num_read_thisfile % BLOCKSIZE == 0))
      Rprintf(" -- Processing %d million reads\n", (num_read_thisfile / BLOCKSIZE + 1) * 10);
    num_read++;
    num_read_thisfile++;

    // only re-process reads withougt perfect hairpin match or without perfect barcode match;
    if ((hairpinindex[num_read] > 0) && (barcodeindex[num_read] > 0))
      continue;

    // re-match barcode:
    if (barcodeindex[num_read] <= 0){
      strncpy(this_barcode, line + barcode_start - 1, barcode_length);
      new_barcode_index = locate_mismatch_barcode(this_barcode);
      if (new_barcode_index > 0)
	barcodecount++;
    } else {
        new_barcode_index = barcodeindex[num_read];
    }
 
    // re-match hairpin:
    if (hairpinindex[num_read] <= 0){
    	strncpy(this_hairpin, line + hairpin_start - 1, hairpin_length);
    	new_hairpin_index = locate_mismatch_hairpin(this_hairpin);
	if ((new_hairpin_index <= 0) && (allow_shifting > 0) && (allow_shifted_mismatch > 0)){
	// Check if given hairpin can be mapped to a shifted location. 
	  int index;
	  // check shifting leftwards
	  for (index = 1; index <= shifting_n_base; index++){
	    strncpy(this_hairpin, line + hairpin_start - 1 - index, hairpin_length);
	    new_hairpin_index = locate_mismatch_hairpin(this_hairpin);
	    if (new_hairpin_index > 0)
	      break;
	  }
	  // check shifting rightwards
	  if (new_hairpin_index <= 0){
	    for (index = 1; index <= shifting_n_base; index++){
	      strncpy(this_hairpin, line + hairpin_start - 1 + index, hairpin_length);
	      new_hairpin_index = locate_mismatch_hairpin(this_hairpin);
	      if (new_hairpin_index > 0)
		break;
	    }
	  }
	}
	if (new_hairpin_index > 0)
	  hairpincount++;
        } else {
    	  new_hairpin_index = hairpinindex[num_read];
        }

      if ((new_barcode_index > 0) && (new_hairpin_index > 0)) {
        summary[new_hairpin_index][new_barcode_index]++;
        bchpcount++;
      }

  }
  fclose(fin);
  free(line);
  free(this_barcode);
  free(this_hairpin);
}


void 
Initialise(int barcodestart, int barcodeend, int hairpinstart, int hairpinend, 
	   int allowshifting, int shiftingnbase,
	   int allowMismatch, int barcodemismatch, int hairpinmismatch, 
	   int allowShiftedMismatch, int verbose){
  int i, j;
  for(i = 0; i < MAX_HAIRPIN; i++) {
    for(j = 0; j < MAX_BARCODE; j++) {
      summary[i][j] = 0;
    }
  }
  num_barcode = 0;
  num_hairpin = 0;

  barcode_start = barcodestart;
  barcode_end = barcodeend;
  hairpin_start = hairpinstart;
  hairpin_end = hairpinend;
  barcode_length = barcode_end - barcode_start + 1;
  hairpin_length = hairpin_end - hairpin_start + 1;

  allow_shifting = allowshifting;
  shifting_n_base = shiftingnbase;
  allow_mismatch = allowMismatch;
  barcode_n_mismatch = barcodemismatch;
  hairpin_n_mismatch = hairpinmismatch;
  allow_shifted_mismatch = allowShiftedMismatch;
  isverbose = verbose;

  num_read = 0;
  barcodecount = 0;
  hairpincount = 0;
  bchpcount = 0;
}

void
Output_Summary_Table(char *output){
  int i, j;
  FILE *fout;
  fout = fopen(output, "w");
  for(i = 1; i <= num_hairpin; i++) {
    fprintf(fout, "%ld", summary[i][1]);
    for(j = 2; j <= num_barcode; j++) {
      fprintf(fout, "\t%ld", summary[i][j]);
    }
    fprintf(fout, "\n");
  }
  fclose(fout);
}

void
Check_Hairpins(void){
  int p, q;
  char base;
  for(p = 1; p <= num_hairpin; p++){
    for(q = 0; q < hairpin_length; q++){
      base = hairpins[p]->sequence[q];
      if ((base != 'A') && (base != 'T') && (base != 'G') && (base != 'C')){
	Rprintf("Hairpin no.%d: %s contains invalid base %c\n", p, hairpins[p]->sequence, base);
      }
    }
  }
}

void
Clean_Up(void){
  int index;
  for (index = 1; index <= num_barcode; index++){
    free(barcodes[index]->sequence);
    free(barcodes[index]);
  }
  for (index = 1; index <= num_hairpin; index++){
    free(hairpins[index]->sequence);
    free(hairpins[index]);
  }
}

void 
processHairpinReads(char **file, int *filecount, 
		    char**barcodeseqs, char**hairpinseqs, 
		    int *barcodestart, int *barcodeend, int *hairpinstart, int *hairpinend, 
		    int *allowShifting, int *shiftingnbase, 
		    int *allowMismatch, int *barcodemismatch, int *hairpinmismatch,
		    int *allowShiftedMismatch,
		    char **output, int *verbose)
{  
  Initialise(*barcodestart, *barcodeend, *hairpinstart, *hairpinend, 
	     *allowShifting, *shiftingnbase,
	     *allowMismatch, *barcodemismatch, *hairpinmismatch, 
	     *allowShiftedMismatch, *verbose);

  Read_In_Barcodes(*barcodeseqs);
  Sort_Barcodes();

  Read_In_Hairpins(*hairpinseqs);
  Check_Hairpins(); 
  Sort_Hairpins();

  long totalreads;
  int i_file;
  totalreads = 0;
  for (i_file = 0; i_file < *filecount; i_file++){
    totalreads = totalreads + Count_Reads(file[i_file]);
  }
  barcodeindex = (int *)malloc(sizeof(int) * totalreads);
  hairpinindex = (int *)malloc(sizeof(int) * totalreads);
  for (i_file = 0; i_file < *filecount; i_file++){
    Process_Hairpin_Reads(file[i_file]);
 }
 
  if (allow_mismatch > 0){
    num_read = 0; 
    // reset total number of read to 0, recheck initial barcode/hairpin index
    Create_Mismatch_Hairpins_List();
    if (num_mismatch_hairpin > 0){
      for (i_file = 0; i_file < *filecount; i_file++){
	Process_Mismatch(file[i_file]);
      }
    }  
  }

  Rprintf("\nThe input run parameters are: \n");
  Rprintf(" -- Barcode: start position %d\t end position %d\t length %d\n", barcode_start, barcode_end, barcode_length);  
  Rprintf(" -- Hairpin: start position %d\t end position %d\t length %d\n", hairpin_start, hairpin_end, hairpin_length); 
  if (allow_shifting) {
    Rprintf(" -- Allow hairpin sequence to be mapped to a shifted position, <= %d base left or right to specified position. \n", shifting_n_base);
  } else {
    Rprintf(" -- Hairpin sequence need to be match at specified position. \n");
  }
  if (allow_mismatch) {
    Rprintf(" -- Allow sequence mismatch, <= %d base in barcode sequence and <= %d base in hairpin sequence. \n", barcode_n_mismatch, hairpin_n_mismatch );
  } else {
    Rprintf(" -- Mismatch in barcode/hairpin sequence is not allowed. \n");
  }

  Rprintf("\nTotal number of read is %ld \n", num_read);
  Rprintf("There are %ld reads (%.4f percent) with barcode match\n", barcodecount, 100.0*barcodecount/num_read);
  Rprintf("There are %ld reads (%.4f percent) with hairpin match\n", hairpincount, 100.0*hairpincount/num_read);
  Rprintf("There are %ld reads (%.4f percent) with both barcode and hairpin match\n", bchpcount, 100.0*bchpcount/num_read);

  Output_Summary_Table(*output);
  free(barcodeindex);
  free(hairpinindex);
}

