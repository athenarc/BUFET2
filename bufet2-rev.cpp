/*
 * Copyright 2021 Konstantinos Zagganas for the Information Management Systems Institute (IMSI) - "Athena" Research Center
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
 * Node used to save the data about each miRNA
 */
struct mNode
{
    string name;
    double left_overlap, center_overlap;
    double left_pvalue, center_pvalue;
    int size;
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
vector<float> i_counts;
vector<long int> left_pvalues,center_pvalues;
go_names_type goNames;
unordered_map<string,string_list> synonyms;
unordered_set<string> genesGo;
vector<mNode *> mList;
/*
 * Global integer variables
 */
int gene_count=1, miRNA_count=1,group_size=0,GOcount=0,thread_count;
bool randomStatus=false;
unsigned long int iterations;
string ontQuery;

/*
 * Function definitions(see explanations above each function)
 */
void getInteractions(string);
void getGOs(string);
void writeOutput(string);
void findPvalues();
void fixInteractions();
void getOntQuery(string);
void calculateCounts();
string trim(string mystr);
string trim_chars_left(string,string);
void getSynonyms(string,string);
bool lcomparison(mNode *, mNode *);
bool ccomparison(mNode *, mNode *);
bool finalComparison(mNode *, mNode *);


int main(int argc, char* argv[])
{   
    
    cout << "Reading GO category data" << endl;
    getGOs(argv[4]);
    if (stoi(argv[7])==0)
    {
        cout << "Reading synonym data" << endl;
        getSynonyms(argv[5],argv[6]);
    }
    else
    {
        cout << "Synonyms disabled" << endl;
    }
    cout << "Reading interaction data" << endl;
    getInteractions(argv[1]);
    
    if (stoi(argv[7])==0)
    {
        cout << "Synonym matching for interactions" << endl;
    }
    fixInteractions();
    cout<< "Reading query term" << endl;
    getOntQuery(argv[3]);
    cout << "Calculating p-values" << endl;
    findPvalues();
    
    cout << "Writing final output" << endl;
    writeOutput(argv[2]);


    return 0;
}

/*
 * The following function reads the differentially expressed miRNAs 
 * in the file provided by the user
 *
 * @param filename: the input file specified by the user
 */
void getOntQuery(string filename)
{
    
    string line,query;
    /*
     * Read file
     */
    ifstream inFile;
    inFile.open(filename);

    if (inFile.is_open())
    {

        getline(inFile,line);
        query=trim(line);
        /*
         * Calculate genes targeted by the given miRNAs
         * 
         * If they do not exist in our hash table
         * increase the counter to be printed as soon as the
         * process ends
         */

        ontQuery=query;
        inFile.close();
    }

    if (ontQuery=="")
    {
        cout << "No ontology terms found in the file. Exiting..." << endl;
    }
    if (goGenes.find(query)==goGenes.end())
    {
        cout << "Ontology term " << query << " was not found in the ontology file" << endl;
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
 * This function calculates the p-values between the query and
 * all miRNAs
 *
 */
void findPvalues()
{
    bits gene_map;
    vec go_map;
            
    go_map=goGenes[ontQuery];
    int total_k=go_map->size();

    for (genedata::iterator it=miRNAs.begin(); it!=miRNAs.end(); it++)
    {   
        gene_map=finalInteractions[it->second];
        int mcount=gene_map->count();
        mNode * node=new mNode;
        node->name=it->first;
        node->size=mcount;
        double intersection=0;
        
        
        
        for (int k=0; k<total_k; k++)
        {
            if ((*gene_map)[(*go_map)[k]]==1)
                    intersection++;
        }
        node->left_overlap=intersection/mcount;
        node->center_overlap=intersection/(mcount+total_k-intersection);

        mList.push_back(node);
            
    }


    double total_m, start, stop;
    /*
     * Get left p-values
     */
    sort(mList.begin(),mList.end(),lcomparison);
    total_m=mList.size();
    stop=total_m-1;
    start=stop-1;

    for(;(stop>=0 && start>=-1);)
    {
        while( (start>=0) && ( (mList[stop]->left_overlap - mList[start]->left_overlap) <1e-7 ) )
        {
            start--;
        }
        /*
         * The for any given overlap the p-value is defined as the number of elements in the table,
         * that are greater or equal to the given sample.
         * This means that if the overlap is the same as some elements in the sorted table, then
         * the p-value is proportional to the index of the element with the lowest index. 
         * So if this index is in position=start+1, the pvalue is total-(position+1)/total
         * 
         * This can never lead to zero p-values.
         */
        double pvalue=(total_m-(start+1))/total_m;
        for (int m=start+1; m<=stop; m++)
        {
            mList[m]->left_pvalue=pvalue;
        }
        stop=start;
        start=stop-1;
        
    }

    /*
     * Get two-sided p-values
     */
    sort(mList.begin(),mList.end(),ccomparison);
    total_m=mList.size();
    stop=total_m-1;
    start=stop-1;

    for(;(stop>=0 && start>=-1);)
    {
        while( (start>=0) && ( (mList[stop]->center_overlap - mList[start]->center_overlap) <1e-7 ) )
        {
            start--;
        }
        /*
         * The for any given overlap the p-value is defined as the number of elements in the table,
         * that are greater or equal to the given sample.
         * This means that if the overlap is the same as some elements in the sorted table, then
         * the p-value is proportional to the index of the element with the lowest index. 
         * So if this index is in position=start+1, the pvalue is total-(position+1)/total
         * 
         * This can never lead to zero p-values.
         */
        double pvalue=(total_m-(start+1))/total_m;
        for (int m=start+1; m<=stop; m++)
        {
            mList[m]->center_pvalue=pvalue;
        }
        stop=start;
        start=stop-1;
        
    }
    
    /*
     * Final sorting using two-sided p-value
     */
    sort(mList.begin(),mList.end(),finalComparison);
}

/*
 * Write final output to a file
 *
 * @param filename: filename provided by the user
 */
void writeOutput(string filename)
{
    ofstream outFile;

    outFile.open(filename);

    if (outFile.is_open())
    {
        /*
         * File header
         */
        outFile << "#miRNA name\tNo-of-miRNA-targets\t";
        outFile << "Left-sided-empirical-p-value\t"; //Benjamini-Hochberg-0.05-FDR" << endl;
        outFile << "Two-sided-empirical-p-value\t";
        outFile << endl;
        for (long unsigned int mi=0; mi < mList.size() ; mi++)
        {
            outFile << mList[mi]->name << "\t" << mList[mi]->size << "\t";
            outFile << mList[mi]->left_pvalue << "\t";
            outFile << mList[mi]->center_pvalue << "\t";
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

bool lcomparison(mNode * n1, mNode * n2)
{
    return (n1->left_overlap < n2->left_overlap);
}
bool ccomparison(mNode * n1, mNode * n2)
{
    return (n1->center_overlap < n2->center_overlap);
}
bool finalComparison(mNode * n1, mNode * n2)
{
    return (n1->center_pvalue < n2->center_pvalue);
}