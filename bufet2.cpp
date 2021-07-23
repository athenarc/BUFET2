/*
 * Copyright 2020 Konstantinos Zagganas for the Information Management Systems Institute (IMSI) - "Athena" Research Center
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * For questions regarding this program, please contact
 * Konstantinos Zagganas at the following e-mail address:
 * zagganas@athenarc.gr
 */

#include <unordered_map>
#include <unordered_set>
#include <string>
#include <iostream>
#include <fstream>
#include <cstring>
#include <random>
#include <vector>
#include <thread>
#include <cstdio>
#include <bitset>
#include <algorithm>

#define BSIZE 25000

using namespace std;

/*
 *
 *
 * Data types defined below
 *
 *
 */

/*
 * Interaction node used for temporary storage of interations
 */
struct interaction
{
    string gene;
    interaction * next;
};

/*
 * GO category node used for temporary storage of GO-gene associations
 */
struct goCatNode
{
    string category;
    string name;
    long int intersection=0;
    long int go_size=0;
    double left_overlap_proportion=0;
    double center_overlap_proportion=0;
    double mean_left_overlap=0;
    double mean_center_overlap=0;
    double left_pvalue=0;
    double center_pvalue=0;
};


/*
 * Hash table containing gene to internal id associations
 */
typedef unordered_map <string,int> genedata;
/*
 * Temporary hash table containing mirna-gene interactions as a linked list for each mirna
 */
typedef unordered_map <int, interaction * > geneInteractions;
/*
 * Temporary hash table containing GO-gene assosiations as a linked list for each category
 */
typedef unordered_map <string, interaction * > goTemp;
/*
 * Auxiliary pointer to a bitset
 */
typedef bitset<BSIZE> * bits;
/*
 * Auxiliary pointer to a vector
 */
typedef vector<int> * vec;
/*
 * Hash table containg interactions for each miRNA
 *
 * @key: internal miRNA id
 * @value: bit vector representing all internal gene ids in set
 */
typedef unordered_map <int, bits> final_interactions_type;
/*
 * Hash table containg GO-gene associations
 *
 * @key: go category
 * @value: bit vector representing all internal gene ids in set
 */
typedef unordered_map <string, vec> goGenes_type;
/*
 * List of tokens
 */
typedef vector<string> token_list;
/*
 * Hash table to contain names for GO categories
 */
typedef unordered_map <string,string> go_names_type;

typedef vector<string> * string_list;

/*
 * Data structures
 */

genedata genes, miRNAs;
geneInteractions interactions;
goGenes_type goGenes;
final_interactions_type finalInteractions;
bits * map_all;
vector<goCatNode> checkGO, noCheckGO;
vector<float> i_counts;
go_names_type goNames;
unordered_map<string,string_list> synonyms;
unordered_set<string> genesGo;

/*
 * Global integer variables
 */
int gene_count=1, miRNA_count=1,group_size=0,GOcount=0,thread_count;
bool randomStatus=false;
unsigned long int iterations;

/*
 * Function definitions(see explanations above each function)
 */
void getInteractions(string);
void getGOs(string);
void getRandom(int,int,int);
void writeOutput(string,int);
void findIntersectionsLeft(int,int);
void findIntersectionsCenter(int,int);
void findIntersections(int,int);
void fixInteractions();
void getMirnas(string);
void calculateCounts();
string trim(string mystr);
string trim_chars_left(string,string);
void getSynonyms(string,string);
bool prepareRandom(int);




int main(int argc, char* argv[])
{   
    thread_count=atoi(argv[6]);
    thread *t= new thread[thread_count];
    iterations=atoi(argv[5]);
    int pmode=atoi(argv[10]);
    map_all=new bits[iterations];
    cout << "Reading GO category data" << endl;
    getGOs(argv[4]);
    if (stoi(argv[9])==0)
    {
        cout << "Reading synonym data" << endl;
        getSynonyms(argv[7],argv[8]);
    }
    else
    {
        cout << "Synonyms disabled" << endl;
    }
    cout << "Reading interaction data" << endl;
    getInteractions(argv[1]);
    
    if (stoi(argv[9])==0)
    {
        cout << "Synonym matching for interactions" << endl;
    }
    fixInteractions();
    cout<< "Calculating query GO overlap" << endl;
    getMirnas(argv[3]);
    cout << "Getting Random miRNA groups" << endl;
    /*
     * Spawn multiple threads to calculate unions
     */
    randomStatus=prepareRandom(group_size);
    if ((thread_count>1) && (!randomStatus))
    {
        if ((group_size==1) && (iterations< finalInteractions.size()))
        {
            getRandom(group_size,0,1);
        }
        else
        {
            for (int u=0; u<thread_count; u++)
                t[u]=thread(getRandom,group_size,u,thread_count);
        
            for (int u=0; u<thread_count; u++)
                t[u].join();
        }
    }
    else
        getRandom(group_size,0,1);
    calculateCounts();
    cout << "Getting GO overlap for " << iterations << " random miRNA groups" << endl;
    /*
     * Spawn multiple threads to calculate intersections
     */
    if (thread_count>1)
    {
        for (int u=0; u<thread_count; u++)
        {
            if (pmode==0)
            {
                t[u]=thread(findIntersections,u,thread_count);
            }
            else if (pmode==1)
            {
                t[u]=thread(findIntersectionsLeft,u,thread_count);
            }
            else if (pmode==2)
            {
                t[u]=thread(findIntersectionsCenter,u,thread_count);
            }
        }
        
        for (int u=0; u<thread_count; u++)
            t[u].join();
    }
    else
    {
        if (pmode==0)
        {
            findIntersections(0,1);
        }
        else if (pmode==1)
        {
            findIntersectionsLeft(0,1);
        }
        else if (pmode==2)
        {
            findIntersectionsCenter(0,1);;
        }
        
    }
    
    cout << "Writing final output" << endl;
    writeOutput(argv[2],pmode);


    return 0;
}

/*
 * The following function reads the differentially expressed miRNAs 
 * in the file provided by the user
 *
 * @param filename: the input file specified by the user
 */
void getMirnas(string filename)
{
    
    string line,miRNA;
    int notFound=0,found=0, target_genes;
    double intersection;
    bitset<BSIZE> gene_map;
    /*
     * Read file
     */
    ifstream inFile;
    inFile.open(filename);

    if (inFile.is_open())
    {

        while(getline(inFile,line))
        {
            miRNA=trim(line);
            /*
             * Calculate genes targeted by the given miRNAs
             * 
             * If they do not exist in our hash table
             * increase the counter to be printed as soon as the
             * process ends
             */
            
            if (miRNA=="") continue;
            if (miRNAs.find(miRNA)!=miRNAs.end())
            {
                found++;
                gene_map |= (*finalInteractions[miRNAs[miRNA]]);

            }
            else
            {
                notFound++;
            }
        }

        inFile.close();
        cout << "Found " << found << " differentially expressed miRNAs" << endl;
        if (notFound)
            cout << notFound << " miRNAs were not found in the set of interactions" << endl;
        group_size=found;
        target_genes=gene_map.count();
        /*
         * Calculate GO categories associated with given miRNAs.
         * Those are the candidates for which to calculate an empirical p-value.
         * The rest have a pvalue=1.
         *
         * Add the results to a special vector
         */
        for (goGenes_type::iterator git=goGenes.begin(); git!=goGenes.end(); git++)
        {   
            /*
             * If the miRNA set has any common genes with a go category, then 
             * add it to the list as a candidate or else add it to a special list
             */
            int total_go=git->second->size();
            intersection=0;
            
            for (int k=0; k<total_go; k++)
            {
                if (gene_map[(*(git->second))[k]]==1)
                    intersection++;
            }
            
            if (intersection>0)
            {
                goCatNode newGO;
                
                newGO.category=git->first;
                newGO.go_size= git->second->size();
                // newGO.intersection=intersection;
                newGO.left_overlap_proportion= intersection/target_genes;
                newGO.center_overlap_proportion= intersection/(total_go+target_genes-intersection);
                newGO.name=goNames[git->first];
                checkGO.push_back(newGO);
                GOcount++;

            }
            else
            {
                goCatNode newGO;

                newGO.category=git->first;
                newGO.go_size= git->second->size();
                newGO.left_overlap_proportion= 0.0;
                newGO.center_overlap_proportion= 0.0;
                newGO.name=goNames[git->first];
                noCheckGO.push_back(newGO);
            }

        }
    }

}


/*
 * This function reads interactions from a file provided by the user.
 *
 * It saves the interactions in a temporary hash table of linked lists
 * @param filename : the file name as provided by the user
 */
void getInteractions(string filename)
{
    
    string line,gline;
    ifstream inFile;
    int index;
    
    inFile.open(filename);
    
    if (inFile.is_open())
    {
        while (getline(inFile,line))
        {
            string miRNA, gene;
            token_list tokens,gtokens;

            line=trim(line);
            if (line=="")
                continue;
            if (line[0]=='#')
                continue;
            
            index=line.find_first_of("|");
            miRNA=line.substr(0,index);
            gene=line.substr(index+1);
            
            /* 
             * If miRNA does not exist in the respective hash tables
             * assign it an internal id and add it
             */         
            
            if (miRNAs.find(miRNA) == miRNAs.end())
            {
                miRNAs[miRNA]=miRNA_count++;
            }
            if (interactions.find(miRNAs[miRNA])==interactions.end())
            {
                interactions[miRNAs[miRNA]]=nullptr;
            }
            
            /*
             * Add interaction to temporary list
             */
            
            interaction * newInteraction= new interaction();

            newInteraction->gene=gene;
            newInteraction->next=interactions[miRNAs[miRNA]];
            interactions[miRNAs[miRNA]]=newInteraction;
        }
        inFile.close();
        
    }
}


/*
 * This function reads GO-gene associations as given by the user.
 *
 * Genes that do not exist in the interactions provided by the user are 
 * assigned an internal id and added to the gene hash table 
 *
 * @param filename: filename provided by the user
 */
void getGOs(string filename)
{
    
    string line;
    ifstream inFile;
    goTemp go_tmp;
    inFile.open(filename);
    unordered_map<string,unordered_set<int> *> tempGO;
    
    /*
     * Temporary saving of data in a hash table containing linked lists
     * like we did for the interactions
     */
    if (inFile.is_open())
    {
        string category, gene, name,domain;
        unsigned int index;

        while (getline(inFile,line))
        {
            index=0;
            
            line=trim(line);
            if (line=="")
                continue;
            if (line[0]=='#')
                continue;
            
            index=line.find_first_of("|");
            gene=line.substr(0,index);
            line=line.substr(index+1);
            index=line.find_first_of("|");
            category=line.substr(0,index);
            line=line.substr(index+1);
            name=line;
            
            

            if (genes.find(gene) == genes.end())
            {
                genes[gene]=gene_count++;
            }
            
            /*
             * Add name to GO names
             */
            goNames[category]=name;
            
            /*
             * If category does not exist add it
             */ 
            if (tempGO.find(category)==tempGO.end())
            {
                tempGO[category]= new unordered_set<int>;
            }
            
            tempGO[category]->insert(genes[gene]);
        }
        inFile.close();
        
        /*
         * Create a vector for each category's genes
         */
        for (unordered_map<string,unordered_set<int> *>::iterator it=tempGO.begin(); it!=tempGO.end(); it++)
        {
            if (goGenes.find(it->first)==goGenes.end())
            {
                    goGenes[it->first]= new vector <int>;
            }
            for (unordered_set<int>::iterator git=(it->second)->begin(); git!=(it->second)->end();git++)
            {   
            
                goGenes[it->first]->push_back((*git));
            }
            delete it->second;
        }
    }
}


/*
 * This function creates a list of synonyms for each gene name
 * 
 * @param filename : filename provided by the user
 */
void getSynonyms(string filename,string taxid)
{
    ifstream inFile;
    string line,gene, taxonomy, synline,synonym;
    int index=0;

    inFile.open(filename);
    if (inFile.is_open())
    {

        while(getline(inFile,line))
        {
            line=trim(line);
            if (line=="")
                continue;
            if (line[0]=='#')
                continue;
            
            index=line.find_first_of("\t");
            taxonomy=line.substr(0,index);
            if (taxonomy!=taxid)
                continue;
            line=line.substr(index+1);
            index=line.find_first_of("\t");
            line=line.substr(index+1);
            index=line.find_first_of("\t");
            gene=line.substr(0,index);
            line=line.substr(index+1);
            index=line.find_first_of("\t");
            line=line.substr(index+1);
            index=line.find_first_of("\t");
            synline=line.substr(0,index);
            
            
            if (synline!="-")
            {
                string_list synonyms_tmp=new vector<string>;

                synonyms_tmp->push_back(trim(gene));
                while((index=synline.find_first_of("|"))!=-1)
                {
                    synonym=synline.substr(0,index);
                    synline=synline.substr(index+1);
                    synonyms_tmp->push_back(trim(synonym));
                }
                /*
                 * Do not forget to add the last synonym
                 */
                synonym=synline;


                synonyms_tmp->push_back(trim(synonym));
            

                for (unsigned int k=0; k<synonyms_tmp->size(); k++)
                {
                    if (synonyms.count((*synonyms_tmp)[k])==0)
                        synonyms[(*synonyms_tmp)[k]]=synonyms_tmp;
            
                }
            }

        }
    }

}

/*
 * This function creates bitsets for each miRNA in the set of interactions and does the basic gene matching
 * of the original empiricalGO.py script
 */
void fixInteractions()
{
    for (geneInteractions::iterator it=interactions.begin(); it!=interactions.end(); it++)
    {

        for (interaction * oldInt=it->second;oldInt!=nullptr; oldInt=oldInt->next)
        {
            string gene=oldInt->gene;
            if (genes.find(gene)==genes.end())
            {
                vector<string> * alternatives;
                if (synonyms.find(gene)!=synonyms.end())
                {
                    alternatives=synonyms[gene];
                }
                else
                {
                    alternatives=nullptr;
                }
                if (alternatives!=nullptr)
                {
                    for (unsigned int k=0; k<alternatives->size(); k++)
                    {
                        if (genes.find((*alternatives)[k])!=genes.end())
                        {
                            oldInt->gene=(*alternatives)[k];
                            break;
                        }
                    }
                }
            }
        }
    }
    for (geneInteractions::iterator it=interactions.begin(); it!=interactions.end(); it++)
    {
        
        finalInteractions[it->first]= new bitset<BSIZE>;

        for (interaction * oldInt=it->second;oldInt!=nullptr; oldInt=oldInt->next)
        {
            string gene=oldInt->gene;
            if (genes.find(gene)==genes.end())
            {
                genes[gene]=gene_count++;
            }
            (*finalInteractions[it->first])[genes[gene]]=1;
        }
    }

}

/* 
 * Check whether the input size is 1.
 * If so the random groups become the 
 * gene sets of the miRNAs in the dataset.
 * 
 * @param size: size of the random miRNA sets
 */
bool prepareRandom(int size)
{
    unsigned long int total_interactions=finalInteractions.size();

    if ((size==1) && (total_interactions < iterations))
    {
        int j=0;
        cout << "You have selected " << iterations << " iterations, but your query contains only 1 miRNA. As a result, only " 
                                     << total_interactions << " iterations can be performed." << endl;
        iterations=total_interactions;
        for (final_interactions_type::iterator fit=finalInteractions.begin(); fit!=finalInteractions.end(); fit++, j++)
        {
            map_all[j]=fit->second;
        }
        return true;

    }
    else
    {
        return false;
    }
}

/*
 * Get a number of random miRNA sets of size m and calculate their interactions 
 * using bitwise operations.
 * The number of random groups is specified by the "iterations" variable
 * The result  will be saved in a vector of bitsets
 * Internal IDs are used
 *
 * @param size: size of the random miRNA sets
 * @param t_num: the thread number
 * @param inc : increment step (the number of threads to be used)
 */

void getRandom(int size,int t_num, int inc)
{
    /*
     * initialize random number generator
     */
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> randMirna(1,miRNA_count-1);

    
    for (unsigned long int i=t_num; i<iterations; i+=inc)
    {
        
        /*
         * Get n random internal IDs where n=size
         * 
         * Bitwise OR for to calculate targeted genes
         */
        bits gene_map=new bitset<BSIZE>;
        
        for (int j=0; j<size; j++)
        {
            int id;
            id=randMirna(gen);
            (*gene_map) |=(*finalInteractions[id]);
        }
        map_all[i]=gene_map;
    }
    
}

/*
 * This function calculates the intsections for all random miRNA sets
 * for the candidate GO categories and calculate the sets with greater overlap
 * than the queried one
 *
 * @param t_num: the thread number
 * @param inc : increment step (the number of threads to be used)
 */
void findIntersections(int t_num,int inc)
{
    bits gene_map;
    vec go_map;
    /*
     * Depending on the thread select specific GO categories to check
     * without overlap between cores
     */
    long int total=checkGO.size();

    for (int git=t_num; git<total; git+=inc)
    {   
        double left_pvalue=0,center_pvalue=0;
        go_map=goGenes[checkGO[git].category];
        int total_k=go_map->size();

        for (unsigned long int i=0; i<iterations; i++)
        {   
            bitset<BSIZE> result;
            double intersection=0;

            gene_map=map_all[i];
            
            for (int k=0; k<total_k; k++)
            {
                if ((*gene_map)[(*go_map)[k]]==1)
                        intersection++;
            }
            double left_overlap=intersection/i_counts[i];
            double center_overlap=intersection/(i_counts[i]+total_k-intersection);

            checkGO[git].mean_left_overlap+=left_overlap;
            checkGO[git].mean_center_overlap+=center_overlap;
            if (left_overlap >= checkGO[git].left_overlap_proportion)
            {   
                left_pvalue++;
            }
            if (center_overlap >= checkGO[git].center_overlap_proportion)
            {   
                center_pvalue++;
            }
            
        }
        checkGO[git].left_pvalue=left_pvalue;
        checkGO[git].center_pvalue=center_pvalue;
    }
}

/*
 * This function calculates the intsections for all random miRNA sets
 * for the candidate GO categories and calculate the sets with greater overlap
 * than the queried one
 *
 * @param t_num: the thread number
 * @param inc : increment step (the number of threads to be used)
 */
void findIntersectionsLeft(int t_num,int inc)
{
    bits gene_map;
    vec go_map;
    /*
     * Depending on the thread select specific GO categories to check
     * without overlap between cores
     */
    long int total=checkGO.size();

    for (int git=t_num; git<total; git+=inc)
    {   
        double left_pvalue=0;
        go_map=goGenes[checkGO[git].category];
        int total_k=go_map->size();

        for (unsigned long int i=0; i<iterations; i++)
        {   
            bitset<BSIZE> result;
            double intersection=0;

            gene_map=map_all[i];
            
            for (int k=0; k<total_k; k++)
            {
                if ((*gene_map)[(*go_map)[k]]==1)
                        intersection++;
            }
            double left_overlap=intersection/i_counts[i];

            checkGO[git].mean_left_overlap+=left_overlap;

            if (left_overlap >= checkGO[git].left_overlap_proportion)
            {   
                left_pvalue++;
            }
            
        }
        checkGO[git].left_pvalue=left_pvalue;
    }
}

/*
 * This function calculates the intsections for all random miRNA sets
 * for the candidate GO categories and calculate the sets with greater overlap
 * than the queried one
 *
 * @param t_num: the thread number
 * @param inc : increment step (the number of threads to be used)
 */
void findIntersectionsCenter(int t_num,int inc)
{
    bits gene_map;
    vec go_map;
    /*
     * Depending on the thread select specific GO categories to check
     * without overlap between cores
     */
    long int total=checkGO.size();

    for (int git=t_num; git<total; git+=inc)
    {   
        double left_pvalue=0,center_pvalue=0;
        go_map=goGenes[checkGO[git].category];
        int total_k=go_map->size();

        for (unsigned long int i=0; i<iterations; i++)
        {   
            bitset<BSIZE> result;
            double intersection=0;

            gene_map=map_all[i];
            
            for (int k=0; k<total_k; k++)
            {
                if ((*gene_map)[(*go_map)[k]]==1)
                        intersection++;
            }
            
            double center_overlap=intersection/(i_counts[i]+total_k-intersection);

            checkGO[git].mean_center_overlap+=center_overlap;
            
            if (center_overlap >= checkGO[git].center_overlap_proportion)
            {   
                center_pvalue++;
            }
            
        }
        checkGO[git].center_pvalue=center_pvalue;
    }
}


/*
 * Write final output to a file
 *
 * @param filename: filename provided by the user
 */
void writeOutput(string filename, int pmode)
{
    ofstream outFile;
    int total_check=checkGO.size(), total_n_check=noCheckGO.size();

    outFile.open(filename);

    if (outFile.is_open())
    {
        /*
         * File header
         */
        outFile << "#GO-term-ID\tGO-term-size\t";
        outFile << "Observed-Target-Left-Sided-Overlap-Proportion\tMean-Random-Simulated-Left-Sided-Overlap-Proportion\tLeft-sided-empirical-p-value\t"; //Benjamini-Hochberg-0.05-FDR" << endl;
        outFile << "Observed-Target-Two-Sided-Overlap-Proportion\tMean-Random-Simulated-Two-Sided-Overlap-Proportion\tTwo-sided-empirical-p-value";
        outFile << endl;
        for (int i=0; i< total_check ; i++)
        {
            outFile << checkGO[i].category << "~" << checkGO[i].name << "\t";
            outFile << checkGO[i].go_size << "\t";
            if (pmode!=2)
            {
                outFile << checkGO[i].left_overlap_proportion << "\t" << checkGO[i].mean_left_overlap/iterations << "\t" << checkGO[i].left_pvalue/iterations << "\t";
            }
            else
            {
                outFile << "-\t-\t-\t";
            }
            if (pmode!=1)
            {
                outFile << checkGO[i].center_overlap_proportion << "\t" << checkGO[i].mean_center_overlap/iterations << "\t" << checkGO[i].center_pvalue/iterations;
            }else
            {
                outFile << "-\t-\t-";
            }

            outFile << endl;
        }
        
        for (int i=0; i < total_n_check ; i++)
        {
            outFile << noCheckGO[i].category << "~" << noCheckGO[i].name << "\t";
            outFile << noCheckGO[i].go_size << "\t";
            if (pmode!=2)
            {
                outFile << "0\t-\t1.0\t";
            }
            else
            {
                outFile << "-\t-\t-\t";
            }
            if (pmode!=1)
            {
                outFile << "0\t-\t1.0\t";
            }
            else
            {
                outFile << "-\t-\t-\t";
            }
            outFile << endl;
        }
        outFile.close();
    }
}

/*
 * This function calculates the number of genes targeted by each random miRNA group
 */
void calculateCounts()
{
    for (unsigned long int i=0; i< iterations; i++)
    {
        i_counts.push_back((float)(map_all[i])->count());
    }
}

/*
 * Auxiliary functions
 */
string trim(string mystr)
{
    int start,stop;

    start=mystr.find_first_not_of("\n\t ");
    stop=mystr.find_last_not_of("\n\t ");

    if ((start==-1) || (stop==-1))
        return "";
    else
        return mystr.substr(start,stop-start+1);
    

}

string trim_chars_left(string mystr,string chars)
{
    int start;

    start=mystr.find_first_not_of(chars);

    return mystr.substr(start);
}