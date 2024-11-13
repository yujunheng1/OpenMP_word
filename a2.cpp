#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <string>
#include <wchar.h>
#include <omp.h>  
#include <memory>
#include <cmath>
#include <locale.h>
#include <fstream>
#include <chrono>
#include <iostream>
#include <locale>
#include <codecvt>
#include <fstream>
#include "hash_table_lock.h"
#define ini_storage 1024
#define byte_letter 4
#define LF 0x0000000A
#define End 0x00000000
#define MAX_WORDS 1000000
#define threshold 10000
#define dict_word 1000000000
#define MAX_CHAR_PER_WORD 1000
#define DICT_SIZE 100000
#define MAX_CANDIDATES 100000
using namespace std;

struct WordRange {
    int start;   // start index
    int end;     // end index
};

int create_length_table(wchar_t** words,int word_size, WordRange length_table[], int max_length) {

    int current_length = wcslen(words[0]);

     for (int i = 0; i <= max_length; ++i) {
        length_table[i].start=-1;  // -1 not length have
        length_table[i].end=-1;
    }

    length_table[current_length].start = 0;

    for (int i = 1; i < word_size; ++i) {
        int word_length = wcslen(words[i]);
        if (word_length != current_length) {
            length_table[current_length].end = i ;

            current_length = word_length;
            length_table[current_length].start = i;
        }
    }

    length_table[current_length].end = word_size ;

    return current_length;
}
int find_start_index(WordRange length_table[],int max_len,int single_word_len){//find the start index based on the length
    for(int i=single_word_len-1;i<=min(single_word_len+1,max_len);i++){
        if(length_table[i].start!=-1){
            return length_table[i].start;
        }
    }
    return -1;
}
int find_end_index(WordRange length_table[],int max_len,int single_word_len){//find the end index based on the length
    for(int i=min(single_word_len+1,max_len);i>=single_word_len-1;i--){
            if(length_table[i].end!=-1){
                return length_table[i].end;
            }
        }
        return -1;
}



void clear(wchar_t** dicts, size_t dicts_size) {
   
    for (size_t i = 0; i < dicts_size; ++i) {
        free(dicts[i]); 
        dicts[i] = nullptr; 
    }
    free(dicts);
    dicts_size=0;
}

int compare_utf32_strings(const wchar_t *str1, const wchar_t *str2) {
    return wcscoll(str1, str2);  
}

void merge_dict(wchar_t **array, size_t left, size_t right) {
    size_t n1 = right - left ; // Length of the left subarray
    wchar_t **tmp = (wchar_t **)malloc((n1 + 1) * sizeof(wchar_t *));

    // Copy data to temporary subarray
    for (size_t i = 0; i <=n1; i++) {
        tmp[i] = array[left + i]; // Use memcpy or assignment
    }

    size_t i = 0, j = (n1/2)+1, k = left;

    // Merge the two halves
   while (i <=(n1/2) && j <=n1)
    {
        // Modify comparison here as needed (e.g., based on key or count)
        int cmp=wcslen(tmp[i])-wcslen(tmp[j]);
        //printf("%ls %ls %zu %zu %d\n",tmp[i].key, tmp[j].key,tmp[i].count,tmp[j].count,cmp);
        if (cmp < 0) {
            array[k] = tmp[i++];
        } else if(cmp>0){
            array[k] = tmp[j++];
        }
        else if(cmp==0){
            int cmp_cnt =compare_utf32_strings(tmp[i], tmp[j]) ;
            if (cmp_cnt < 0) {
                array[k] = tmp[i++];
            } else{
                array[k] = tmp[j++];
            }
        }
         k++;
    }
     while (i <=n1/2){
        array[k++] = tmp[i++];
    }

    // Copy any remaining elements of rightarr
    while (j <=n1){
        array[k++] = tmp[j++];
    }
    free(tmp);  
}

void parallel_merge_sort_dict(wchar_t **array, size_t left, size_t right) {
    if (left < right) {
        size_t mid = left + (right - left) / 2;

        if ((right - left) > threshold) {
            #pragma omp task shared(array)
            parallel_merge_sort_dict(array, left, mid);

            #pragma omp task shared(array)
            parallel_merge_sort_dict(array, mid + 1, right);
            
            #pragma omp taskwait
        } else {
            parallel_merge_sort_dict(array, left, mid);
            parallel_merge_sort_dict(array, mid + 1, right);
        }
        merge_dict(array, left, right);
    }
}


void merge_keyCount(KeyCountPair *array, size_t left, size_t right) {
    size_t n1 = right - left ; // Length of the left subarray
    KeyCountPair *tmp = (KeyCountPair *)malloc((n1+1) * sizeof(KeyCountPair));

    // Copy data to temporary subarray
    for (size_t i = 0; i <=n1; i++) {
        tmp[i] = array[left + i]; // Use memcpy or assignment
    }

    size_t i = 0, j = (n1/2)+1, k = left;

    // Merge the two halves
   while (i <=(n1/2) && j <=n1)
    {
        // Modify comparison here as needed (e.g., based on key or count)
        int cmp = tmp[i].count-tmp[j].count;
        //printf("%ls %ls %zu %zu %d\n",tmp[i].key, tmp[j].key,tmp[i].count,tmp[j].count,cmp);
        if (cmp < 0) {
            array[k] = tmp[i++];
        } else if(cmp>0){
            array[k] = tmp[j++];
        }
        else if(cmp==0){
            int cmp_cnt =compare_utf32_strings(tmp[i].key, tmp[j].key) ;
            if (cmp_cnt < 0) {
                array[k] = tmp[i++];
            } else{
                array[k] = tmp[j++];
            }
        }
         k++;
    }
     while (i <=n1/2){
        array[k++] = tmp[i++];
    }

    // Copy any remaining elements of rightarr
    while (j <=n1){
        array[k++] = tmp[j++];
    }
    free(tmp);  
}

void parallel_merge_sort_keyCount(KeyCountPair *array, size_t left, size_t right) {
    if (left < right) {
        size_t mid = left + (right - left) / 2;

        if ((right - left) > threshold) {
            #pragma omp task shared(array)
            parallel_merge_sort_keyCount(array, left, mid);

            #pragma omp task shared(array)
            parallel_merge_sort_keyCount(array, mid + 1, right);
            
            #pragma omp taskwait
        } else {
            parallel_merge_sort_keyCount(array, left, mid);
            parallel_merge_sort_keyCount(array, mid + 1, right);
        }
        merge_keyCount(array, left, right);
    }
}

int dict_chunk_proc(int size,int rank, const char* filename,size_t file_size,wchar_t*** wordList,size_t chunk_size,int &start_chunk,int &End_chunk)
 {//use for the first time to access the file and need to find the start and end of the word; show broadcasts for each node
    
    size_t start = rank * chunk_size*byte_letter;
    size_t end_chunk = (rank == size - 1) ? (file_size-byte_letter) : (rank + 1) * chunk_size * byte_letter;
    FILE *file;


    file = fopen(filename, "rb");

    fseek(file, start, SEEK_SET);

    if (start != 0) {//find the LF the start of a new word
        wchar_t c;
        
        fseek(file, start, SEEK_SET); 

        while (start > 0) {
            c = fgetc(file) | (fgetc(file) << 8) | (fgetc(file) << 16) | (fgetc(file) << 24);
            if (c == LF) { 
                (start) += byte_letter;  
                fseek(file, start, SEEK_SET);
                break;
            }
            
            (start) += byte_letter;  
            if(start>=(file_size-byte_letter)){
                *wordList=NULL;
                return 0;
            }
            fseek(file, start, SEEK_SET); 
        }
    }
    else 
    {
        start=start+4;
        fseek(file, start, SEEK_SET); 
    }
    
    start_chunk=start;
    wchar_t **words = (wchar_t **)malloc(MAX_WORDS * sizeof(wchar_t *));
    
    wchar_t temp_word[MAX_CHAR_PER_WORD + 1];  
    
    wchar_t c;
    int index=0;
    size_t char_index=0;
    while (1)
    {
        c = fgetc(file) | (fgetc(file) << 8) | (fgetc(file) << 16) | (fgetc(file) << 24);
        if ((start >= end_chunk&&c == LF)) {//finish the end word of current chunk
            temp_word[char_index] = End;
            words[index] = (wchar_t *)malloc((char_index + 1) * sizeof(wchar_t));
            wcscpy(words[index], temp_word);
            index++;
            break;
        }
        if (c == LF) {
            if (char_index > 0) {
                temp_word[char_index] = End;
                words[index] = (wchar_t *)malloc((char_index + 1) * sizeof(wchar_t));
                wcscpy(words[index], temp_word);
                index++;
                char_index = 0; 
            }
            (start) += byte_letter;
            continue; 
        }
        temp_word[char_index] = c;
        char_index++; 
        (start) += byte_letter;
        fseek(file, start, SEEK_SET);
    }
    End_chunk=start;
    *wordList=words;
    return index;

 }

 int dict_chunk_proc_have_range(int size,int rank, const char* filename,wchar_t*** wordList,int start_chunk,int end_chunk)
 {//processing the file with corrent index without search the start of word and the end of the word
    FILE *file;
    size_t start=start_chunk;
    file = fopen(filename, "rb");

    fseek(file, start_chunk, SEEK_SET);

    
    wchar_t **words = (wchar_t **)malloc(MAX_WORDS * sizeof(wchar_t *));
    
    wchar_t temp_word[MAX_CHAR_PER_WORD + 1];  
    
    wchar_t c;
    int index=0;
    size_t char_index=0;
    while (1)
    {
        c = fgetc(file) | (fgetc(file) << 8) | (fgetc(file) << 16) | (fgetc(file) << 24);
        
        if (start == static_cast<size_t>(end_chunk)) {
          //  printf("bre:%zu %zu\n",index,char_index);
            temp_word[char_index] = End;
            words[index] = (wchar_t *)malloc((char_index + 1) * sizeof(wchar_t));
            wcscpy(words[index], temp_word);
            index++;
            break;
        }
        if (c == LF) {
            if (char_index > 0) {
                temp_word[char_index] = End;
                words[index] = (wchar_t *)malloc((char_index + 1) * sizeof(wchar_t));
                wcscpy(words[index], temp_word);
                index++;
                char_index = 0; 
            }
            (start) += byte_letter;
            continue; 
        }
        temp_word[char_index] = c;
        char_index++;
        (start) += byte_letter;
        fseek(file, start, SEEK_SET);
    }
    end_chunk=start;
   
    *wordList=words;
    
    return index;

 }


int split_chunk(int size,int rank, const char* filename,wchar_t*** wordList)
{//use for word file
    
    FILE *file;
    size_t file_size;
    size_t chunk_size, start, end_chunk;
    if (rank == 0) {
        file = fopen(filename, "rb");
        fseek(file, 0, SEEK_END);
        file_size = ftell(file);
        rewind(file);
        
        chunk_size = (file_size/byte_letter) / size;

        fclose(file);
    }
    MPI_Bcast(&file_size, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&chunk_size, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    
    start = rank * chunk_size*byte_letter;
    end_chunk = (rank == size - 1) ? (file_size-byte_letter) : (rank + 1) * chunk_size * byte_letter;

\
    file = fopen(filename, "rb");

    fseek(file, start, SEEK_SET);

    if (start != 0) {
        wchar_t c;
        
        fseek(file, start, SEEK_SET); 

        while (start > 0) {
            c = fgetc(file) | (fgetc(file) << 8) | (fgetc(file) << 16) | (fgetc(file) << 24);
            if (c == LF) { 
                (start) += byte_letter;  
                fseek(file, start, SEEK_SET);
                break;
            }
            
            (start) += byte_letter;  
            if(start>=(file_size-byte_letter)){
                *wordList=NULL;
                return 0;
            }
            fseek(file, start, SEEK_SET); 
        }
    }
    else 
    {
        start=start+4;
        fseek(file, start, SEEK_SET); 
    }
    
    
    wchar_t **words = (wchar_t **)malloc(MAX_WORDS * sizeof(wchar_t *));
    
    wchar_t temp_word[MAX_CHAR_PER_WORD + 1];  
    
    wchar_t c;
    int index=0;
    size_t char_index=0;
   
    while (1)
    {
        c = fgetc(file) | (fgetc(file) << 8) | (fgetc(file) << 16) | (fgetc(file) << 24);
    
        if ((start >= end_chunk&&c == LF)) {
          
            temp_word[char_index] = End;
            words[index] = (wchar_t *)malloc((char_index + 1) * sizeof(wchar_t));
            wcscpy(words[index], temp_word);
            index++;
            break;
        }
        if (c == LF) {
            if (char_index > 0) {
                temp_word[char_index] = End;
                words[index] = (wchar_t *)malloc((char_index + 1) * sizeof(wchar_t));
                wcscpy(words[index], temp_word);
                index++;
                char_index = 0; 
            }
            (start) += byte_letter;
            continue; 
        }
        temp_word[char_index++] = c;
        (start) += byte_letter;
        fseek(file, start, SEEK_SET);
    }
   
    *wordList=words;
    return index;
    
}

int is_equal(const wchar_t  letterWord,const wchar_t letterDict,bool is_all_caps){
    if(is_all_caps){
        
        return toupper(letterDict)==letterWord;
    }
    return (letterDict==letterWord);
}
int is_equal_cp(const wchar_t shorter,const wchar_t longer){

    
    return (longer==shorter);
    
}

int differsByOneChange(const wchar_t* word, const wchar_t* dict){
    int is_all_caps = 1;
    int i_word=0;

    int diffCount = 0;
    int lenWord = wcslen(word);
    int lenDict = wcslen(dict);
    
    if(abs(lenDict-lenWord)>1){//when the length is large than 1 do not cmp
        return -1;
    }

    for (i_word = 0; word[i_word] != End; i_word++) {
        if (!isupper(word[i_word])) {
            is_all_caps = 0;
            break;
        }
    }
    
    if (lenDict == lenWord) {
        for (int i = 0; i < lenWord; i++) {
            if (is_equal(word[i],dict[i],is_all_caps)==1){
                continue;
            }   
            diffCount++;
                
            if (i < lenWord - 1 && word[i] == dict[i + 1] && word[i + 1] == dict[i]) {
                i++;
            }
                // If more than one letter differs, return false
            if (diffCount > 1) {
                return -1;
            }
        }
        if(diffCount==1&&iswupper(word[0])){
            if (is_equal(tolower(word[0]),dict[0],is_all_caps)==1){
                return 0;
            }
        }
        return diffCount;
    }
    if (abs(lenDict-lenWord) == 1) {
        const wchar_t* longer = lenWord > lenDict  ? word:dict;
        const wchar_t* shorter = lenWord < lenDict ? word:dict;
        int lenLonger = max(lenDict, lenWord);
        int lenShorter = min(lenDict, lenWord);

        int i_short = 0, i_long = 0;
        while (i_short < lenShorter && i_long < lenLonger) {
    
            if(is_equal_cp(shorter[i_short],longer[i_long])==0){
                diffCount++;
                i_long++;  
                if (diffCount > 1) {
                    return -1;  
                }
            }
            else{
                i_long++;
                i_short++;
            }
           
        }
        
        return diffCount+(abs((lenWord+lenDict)-(i_short+i_long)));  // Only one change (insertion or deletion)        }
        
    }
    return diffCount;
}
void freeDictionary(wchar_t** dictionary, size_t size) {
    for (size_t i = 0; i < size; i++) {
        free(dictionary[i]); 
    }
    delete[] dictionary; 
}


int check_if_a_word_in_dictionary(int start_dict, int end_dict,wchar_t* word, wchar_t** dictionary, size_t dict_size,HashTable* hashTable_word) {
    
    int candidateCount = 0;

    wchar_t wordCopy[MAX_CHAR_PER_WORD];
    wchar_t** condiction_word = new wchar_t*[MAX_CANDIDATES];
    
    wcscpy(wordCopy, word);
    for(int i=start_dict;i<end_dict;i++ ){
        wchar_t dictWordCopy[MAX_CANDIDATES];
        wcscpy(dictWordCopy, dictionary[i]);
        int is_add=differsByOneChange(wordCopy, dictWordCopy);
        if(is_add==0){
            freeDictionary(condiction_word,candidateCount);
            return 2;
        }
        if ((is_add==1)) {
            condiction_word[candidateCount]=wcsdup(dictWordCopy);
            candidateCount++;
        }
        if(candidateCount>MAX_CANDIDATES){
            printf("no moemory\n");
        }
    }
    if(candidateCount==0){
        return 0;
    }
    
    for(int i=0;i<candidateCount;i++){
        hashTable_word->insert(word,wcsdup(condiction_word[i]));
    }
    freeDictionary(condiction_word,candidateCount);
    return 1;
}

void deal_word_dict(wchar_t **words, wchar_t **dict,size_t word_cnt, size_t dict_cnt,HashTable* hashTable_word ,int* queried){
    
    int num_thread=max(min(8,(int)log(dict_cnt)),1);
    //merge sort
    #pragma omp parallel num_threads(num_thread)
    {
        #pragma omp single
        {
            parallel_merge_sort_dict(dict,0,dict_cnt-1);
        }
    }
   
    const int max_length = wcslen(dict[dict_cnt-1]);
    WordRange length_table[max_length+1];
    int table_size = create_length_table(dict, dict_cnt,length_table, max_length);

    #pragma omp parallel
    {
        #pragma omp for 
        for(size_t i=0;i<word_cnt;i++){
            
            if(queried[i]==2){// 0 not find; 1 add into condition; 2 in dict
                continue;
            }
            int single_word_len=wcslen(words[i]);
            if(single_word_len-1>table_size){
                continue;
            }
            int start_dict=find_start_index(length_table,table_size,single_word_len);
            int end_dict=find_end_index(length_table,table_size,single_word_len);
            if(start_dict==-1||end_dict==-1){
                continue;
            }
            int is_put=check_if_a_word_in_dictionary(start_dict,end_dict,words[i],dict,dict_cnt,hashTable_word);
        
            if(queried[i]==1&&is_put==2){
                hashTable_word->remove(words[i]);
                queried[i]=2;
            }
            else{//Flags are prioritized base on value
                queried[i]=max(is_put,queried[i]);
            }
            
            
        }
    }
    #pragma omp barrier
}

void getallhash(int rank,int size,HashTable* hashTable_word,MPI_Comm comm){
    if(rank==0){
        char** buffer_list = new char*[size-1];
        size_t* total_size_list = new size_t[size - 1];
        MPI_Request requests1[size-1];  //request for total size
        
        MPI_Request requests2[size-1];  //request for buffer
         
       
        for(int i=1;i<size;i++){
            MPI_Irecv(&total_size_list[i-1], 1, MPI_UNSIGNED_LONG, i, 0, comm, &requests1[i-1]);
            MPI_Wait(&requests1[i-1], MPI_STATUS_IGNORE); 
            char* buffer = new char[total_size_list[i-1]];
            MPI_Irecv(buffer, total_size_list[i-1] , MPI_BYTE, i, 1, comm, &requests2[i-1]);
            MPI_Wait(&requests2[i-1], MPI_STATUS_IGNORE); 
            buffer_list[i-1]=buffer;
        }
        
        #pragma omp parallel num_threads(size-1)
        {
            #pragma omp for
            for(int i=1;i<size;i++)
            {
                
                rec_deserialize(*hashTable_word, buffer_list[i-1],total_size_list[i-1]);
                
            }
        }
    }
    else{
        send_serialize(*hashTable_word,rank,0,comm);
    }
}

void printHashTable(HashTable& hashTable,const KeyCountPair* pairs,const size_t keyCnt) {
    for (size_t i = 0; i < keyCnt; i++) {
        size_t cntValue;
        printf("%ls: ",pairs[i].key);
        wchar_t** value = hashTable.getAllValue(pairs[i].key,cntValue); // Get value from hash table
        for(size_t i=0;i<cntValue;i++){
            printf("%ls,",value[i]);
        }
        printf("\n");
    }
}
void proce_dict_chunk(int rank, int size,const char* filename,wchar_t **words,size_t words_size, HashTable* hashTable_word,int* queried,MPI_Comm comm){
    FILE *file;
    int file_size=0;
    size_t chunk_size;
    int* gobal_chunk_range = new int[2 * size];
    for (int i = 0; i < 2 * size; i++) {
        gobal_chunk_range[i] = -1;
    }
    if(rank==0)
    {
        file = fopen(filename, "rb");
        fseek(file, 0, SEEK_END);
        file_size = ftell(file);
        rewind(file);
        
        chunk_size = (file_size/byte_letter) / size;

        fclose(file);
    }
        MPI_Bcast(&file_size, 1, MPI_LONG, 0, MPI_COMM_WORLD);
        MPI_Bcast(&chunk_size, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    
    int local_start;
    int local_end;

    int dicts_size;
    wchar_t** dicts = NULL; 
    dicts_size=dict_chunk_proc(size,rank,filename,file_size,&dicts,chunk_size,local_start,local_end);
    
    deal_word_dict(words, dicts,words_size, dicts_size,hashTable_word,queried);
    if(rank==0){
        printf("Round 0 Rank 0 have %d words in dictionary \n",dicts_size);
    }
    
    

    int local_data[2];
    local_data[0] = local_start; local_data[1] = local_end;
    MPI_Allgather(local_data, 2, MPI_INT, gobal_chunk_range, 2, MPI_INT, comm);
    
    for(int i=1;i<size;i++){
        int index_chunk=(rank+i)%size;
        int start_chunk=gobal_chunk_range[index_chunk*2];
        int end_chunk=gobal_chunk_range[index_chunk*2+1];
        
        dicts_size=dict_chunk_proc_have_range(size,rank,filename,&dicts,start_chunk,end_chunk);
        if(rank==0){
        printf("Round %d Rank 0 have %d words in dictionary \n",i,dicts_size);
        }
        deal_word_dict(words, dicts,words_size, dicts_size,hashTable_word,queried);
        
        clear(dicts,dicts_size);
        #pragma omp barrier
        
    }

    
}

void writeHashTableToFile(HashTable& hashTable,  KeyCountPair* pairs,  size_t keyCnt,  const char* fileName) {
    FILE* outfile = fopen(fileName, "w");
    fputc(L'\xff', outfile);
    fputc(L'\xfe', outfile);
    fputc(L'\0', outfile);
    fputc(L'\0', outfile);

   

    for (size_t i = 0; i < keyCnt; i++) {
        size_t cntValue;
        
    
        wchar_t** value = hashTable.getAllValue(pairs[i].key, cntValue);
        
       
         wchar_t* key = pairs[i].key;
        std::wstring words_str;
        for (size_t j = 0; j < cntValue; j++) {
           
            words_str += value[j];     
            words_str += L" ";         
        }

        words_str += L"\n";
        std::wstring line_str;
        line_str =  std::wstring(key) + L":" + L" " + words_str;
        fwrite(line_str.c_str(), sizeof(wchar_t), line_str.size(), outfile);
    }

    fclose(outfile);
}


void print_running_time(string str,int rank,std::chrono::high_resolution_clock::time_point start,std::chrono::high_resolution_clock::time_point end){
     // Show wall time
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
  // Convert duration to minutes, seconds, and milliseconds
  auto minutes = duration / 60000;
  auto seconds = (duration % 60000) / 1000;
  auto milliseconds = duration % 1000;
  if(rank!=-1){
  cout << "[*] Wall time: "<<str <<rank<<" "<< minutes << " mins " << seconds << " secs "
       << milliseconds << " ms" << "; " << duration << endl;
    }
    else{
        cout << "[*] All Wall time: "<<str <<" "<< minutes << " mins " << seconds << " secs "
       << milliseconds << " ms" << "; " << duration << endl;
    }
}



int main(int argc, char** argv) {
    
    auto start = std::chrono::high_resolution_clock::now(); 
    MPI_Init(&argc, &argv);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (argc != 3) {
        if (rank == 0) {
            
            printf("Usage: ./spellcheck dictionary.txt words.txt");
        }
        MPI_Finalize();
        return 1;
    }

    
    const char* dictionary_file = argv[1];  
    const char* words_file = argv[2];        

    wchar_t **  words = NULL;

    int words_size;

    HashTable ht(DICT_SIZE);
    
    words_size=split_chunk(size, rank,words_file,&words);
    printf("Process %d have %d words\n",rank,words_size);
   
    int queried[words_size] = {0};

    auto end_word = std::chrono::high_resolution_clock::now(); 
    print_running_time("Word file read ",rank,start,end_word);

    proce_dict_chunk(rank,size,dictionary_file,words,words_size,&ht,queried,MPI_COMM_WORLD);


    

   for(int i=0;i<words_size;i++)
   {
    //printf("new %ls %d\n",words[i],queried[i]);
        if(queried[i]==0){
            ht.insert(wcsdup(words[i]),L"");
        }
    } 

    auto end_word_pro = std::chrono::high_resolution_clock::now(); 
    print_running_time("Word compare ",rank,end_word,end_word_pro);


    size_t keyCount;
    KeyCountPair* pairs = ht.getAllKeysAndCounts(keyCount);
    #pragma omp parallel num_threads(8)
    {
        #pragma omp single
        {
            parallel_merge_sort_keyCount(pairs,0,keyCount-1);
        }
    }

    auto start_tran = std::chrono::high_resolution_clock::now(); 

    getallhash(rank,size,&ht,MPI_COMM_WORLD); 
    MPI_Barrier(MPI_COMM_WORLD);

    auto end_tran = std::chrono::high_resolution_clock::now(); 
    print_running_time("Transfer ",rank,start_tran,end_tran);

   
    
    if(rank==0){
       
        size_t keyCount;
        KeyCountPair* pairs = ht.getAllKeysAndCounts(keyCount);
        #pragma omp parallel num_threads(8)
        {
            #pragma omp single
            {
                parallel_merge_sort_keyCount(pairs,0,keyCount-1);
            }
        }
        //printHashTable(ht,pairs,keyCount);
        std::string fileName = "output.txt";
        writeHashTableToFile(ht,pairs,keyCount,fileName.c_str());

    }
    

    MPI_Finalize();
    
    auto end = std::chrono::high_resolution_clock::now(); 
    print_running_time("Whole Program ",-1,start,end);
    return 0;
}
