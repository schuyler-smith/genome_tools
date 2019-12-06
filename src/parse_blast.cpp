
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iostream>
#include <filesystem>

int main()
{
	std::string 	subject,
					e_value;

	int 			length,
					mismatch,
					gap,
					qstart,
					qend,
					sstart,
					send,
					count;

	double 			perc_id,
					perc_id_cutoff,
					bitscore;

	perc_id_cutoff = 100;
	count = 0;

	std::ifstream 	file("/home/schuyler/Dropbox/read_blast/blast_results/test.blast");

	// std::ifstream 	file("/home/schuyler/Desktop/read_blast/Q3-Mock-0-009spike-A_S7_L001_R1_001.blast");
 	std::string 	line;

    while(std::getline(file, line)) 
    {
    	std::stringstream   line_s(line);
    	std::string 		query;

    	std::getline(line_s, query, '\t');
    	line_s >> subject >> perc_id >> length >> mismatch >> gap >> 
    		qstart >> qend >> sstart >> send >> e_value >> bitscore;
    	if(perc_id >= perc_id_cutoff)
    	{
    		++count;
    		std::cout << "" << query << '\n' ;
    	}
    }

return(count);
}