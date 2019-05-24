package data_create;

/**
 * @author Markus Hoffmann
 * @contact Markus.Hoffmann@campus.lmu.de
 * @input miRNA-seq, exp-seq, miRNA-mapping, cutoff(for donors - optional, if not set 20 is used - on place 5 of args)
 * @output directory for miRNA matrix, gene matrix
 * @information: only miRNAs, genes are used which are found in at least in cutoffs samples
 * @datelastmodified 2019/05/24
 */


import java.io.BufferedReader;
import java.io.FileReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import javax.xml.crypto.Data;

public class DataExtraction {
	
	private static class DataStruc {
	public HashMap<String,Double> genes = new HashMap<>();
	public HashMap<String,Double> miRNAs = new HashMap<>();
	
	public DataStruc()
	{
		
	}
	
	public void addGene(String gene, double expression)
	{
		if(!genes.containsKey(gene))
		{
			genes.put(gene,expression);
		}
	}
	public void addmiRNA(String miRNA, double expression)
	{
		if(!miRNAs.containsKey(miRNA))
		{
			miRNAs.put(miRNA,expression);
		}
	}
}
	
	public static void main(String[] args)
	{
		
		String input_miRNAseq = args[0];
		String input_expSeq = args[1];
		String input_miRNA_map = args[2];
		String output = args[3];
		
		int cutoff= 20;
		if(args.length==5)
		{
			cutoff = Integer.parseInt(args[4]);
		}

		
		try 
		{
			
			HashMap<String, DataStruc> sample_genes_miRNAs = new HashMap<>();
			//unique genes and counts of donors
			HashMap<String,Integer> unique_genes = new HashMap<>();
			//unique miRNAs and counts of donors
			HashMap<String,Integer> unique_miRNAs = new HashMap<>();
			//samples which are already counted for unique miRNA and unique genes -> is cleared after use
			HashSet<String> unique_samples = new HashSet<>();
			//mapping of miRNA hsa codes to MIMAT codes
			HashMap<String,String> miRNA_mapper = new HashMap<>();
			/**
			 * mapping of ENSG numbers to gene symbols could be added easily
			 */
			
			//read and save mapping for miRNA (hsa codes to MIMAT codes)
			BufferedReader br_mapper_miRNA = new BufferedReader(new FileReader(input_miRNA_map));
			String lineFile_miRNAmap ="";
			while((lineFile_miRNAmap = br_mapper_miRNA.readLine()) != null)
			{
				String[] inf = lineFile_miRNAmap.split("\t");
				miRNA_mapper.put(inf[0], inf[1]);
				
			}
			
			br_mapper_miRNA.close();
			
			//read and save miRNA and expression lvls
			BufferedReader br_miRNA = new BufferedReader(new FileReader(input_miRNAseq));
			
			String lineFile = "";
			int col_miRNAexp = 0;
			int col_miRNA = 0;
			while((lineFile = br_miRNA.readLine()) != null)
			{
				String[] inf = lineFile.split("\t");
				if(lineFile.startsWith("icgc"))
				{
					//search for mirna and expression level in header
					for(int i = 0; i < inf.length; i++)
					{
						if(inf[i].equals("mirna_id"))
						{
							col_miRNA = i;
						}
						if(inf[i].equals("normalized_read_count"))
						{
							col_miRNAexp = i;
						}
					}
					continue;
				}
				//map the codes
				String map = miRNA_mapper.get(inf[col_miRNA]);
				//count this miRNA for every donor
				if(unique_miRNAs.containsKey(map)&& !unique_samples.contains(inf[0]+"_"+inf[col_miRNA]))
				{
					int z = unique_miRNAs.get(map);
					z++;
					unique_miRNAs.put(map,z);
				}
				if(!unique_miRNAs.containsKey(map)&& !unique_samples.contains(inf[0]+"_"+inf[col_miRNA]))
				{
					unique_miRNAs.put(map,1);
				}
				unique_samples.add(inf[0]+"_"+inf[col_miRNA]);
				
				//add to data structure of donor
				if(!sample_genes_miRNAs.containsKey(inf[0]))
				{
					DataStruc s = new DataStruc();
					s.addmiRNA(map,Double.parseDouble(inf[col_miRNAexp]));
					
					sample_genes_miRNAs.put(inf[0], s);
				}
				else
				{
					DataStruc s = sample_genes_miRNAs.get(inf[0]);
					s.addmiRNA(map,Double.parseDouble(inf[col_miRNAexp]));
				}
			}
			
			br_miRNA.close();
			unique_samples.clear();
			
			
			//read and save genes with expression levels
			BufferedReader br_expseq = new BufferedReader(new FileReader(input_expSeq));

			String lineExp = "";
			int col_gene = 7;
			int col_genExp =8;
			while((lineExp = br_expseq.readLine()) != null)
			{
				String inf[] = lineExp.split("\t");
				//searching for gene id and expression level in header
				if(lineExp.startsWith("icgc"))
				{
					for(int i = 0; i < inf.length; i++)
					{
						if(inf[i].equals("gene_id"))
						{
							col_gene = i;
						}
						if(inf[i].equals("normalized_read_count"))
						{
							col_genExp=i;
						}
					}
					
					continue;
				}
				//count genes for every donor
				if(unique_genes.containsKey(inf[col_gene])&&!unique_samples.contains(inf[0]+"_"+inf[col_gene]))
				{
					int z = unique_genes.get(inf[col_gene]);
					z++;
					unique_genes.put(inf[col_gene], z);
				}
				if(!unique_genes.containsKey(inf[col_gene])&&!unique_samples.contains(inf[0]+"_"+inf[col_gene]))
				{
					unique_genes.put(inf[col_gene], 1);
				}
				unique_samples.add(inf[0]+"_"+inf[col_gene]);
				
				if(!sample_genes_miRNAs.containsKey(inf[0]))
				{
					DataStruc s = new DataStruc();
					s.addGene(inf[col_gene], Double.parseDouble(inf[col_genExp]));
					
					sample_genes_miRNAs.put(inf[0], s);
				}
				else
				{
					DataStruc s = sample_genes_miRNAs.get(inf[0]);
					s.addGene(inf[col_gene],Double.parseDouble(inf[col_genExp]) );
				}
			}
			
			br_expseq.close();

			//print genes matrix file
			PrintWriter pw_genes = new PrintWriter(output+"\\"+"genes.tsv");

			ArrayList<String> header_al_genes = new ArrayList<>();
			for(String k: unique_genes.keySet())
			{
				
				if(unique_genes.get(k)>=cutoff)
				{
					header_al_genes.add(k);
				}
			}
			
			StringBuilder sb_genes = new StringBuilder();
			sb_genes.append(" ");
			for(int x = 0; x < header_al_genes.size();x++)
			{
				sb_genes.append("\t"+header_al_genes.get(x));
			}
			
			pw_genes.println(sb_genes.toString());
			//System.out.println(sb_genes.toString());
			
			for(String k :sample_genes_miRNAs.keySet())
			{
				StringBuilder sb_lines_gene = new StringBuilder();
				sb_lines_gene.append(k);
				
				DataStruc s = sample_genes_miRNAs.get(k);
				HashMap<String,Double> hm = s.genes;
				for(int x = 0; x < header_al_genes.size(); x++)
				{
					sb_lines_gene.append("\t");
					if(hm.containsKey(header_al_genes.get(x)))
					{
						sb_lines_gene.append(hm.get(header_al_genes.get(x)));
					}
					else
					{
						sb_lines_gene.append("0");
					}
					
				}
				pw_genes.println(sb_lines_gene.toString());
				//System.out.println(sb_lines_gene.toString());
			}
			pw_genes.close();
			
			PrintWriter pw_miRNA = new PrintWriter(output+"\\"+"miRNA.tsv");
			//String[] header_arr_miRNA = new String[unique_miRNAs.size()];
			//unique_miRNAs.toArray(header_arr_miRNA);
			
			ArrayList<String> header_al_miRNA = new ArrayList<>();
			for(String k: unique_miRNAs.keySet())
			{
				if(unique_miRNAs.get(k)>=cutoff)
				{
					header_al_miRNA.add(k);
				}
			}
			
			StringBuilder sb_miRNA = new StringBuilder();
			sb_miRNA.append(" ");
			for(int x = 0; x < header_al_miRNA.size();x++)
			{
				sb_miRNA.append("\t"+header_al_miRNA.get(x));
			}
			
			pw_miRNA.println(sb_miRNA.toString());
			
			for(String k :sample_genes_miRNAs.keySet())
			{
				StringBuilder sb_lines_miRNA = new StringBuilder();
				sb_lines_miRNA.append(k);
				
				DataStruc s = sample_genes_miRNAs.get(k);
				HashMap<String,Double> hm = s.miRNAs;
				for(int x = 0; x < header_al_miRNA.size(); x++)
				{
					sb_lines_miRNA.append("\t");
					if(hm.containsKey(header_al_miRNA.get(x)))
					{
						sb_lines_miRNA.append(hm.get(header_al_miRNA.get(x)));
					}
					else
					{
						sb_lines_miRNA.append("0");
					}
				}
				pw_miRNA.println(sb_lines_miRNA.toString());
			}
			pw_miRNA.close();
		}
		catch(Exception e)
		{
			e.printStackTrace();
		}
		
	}
		

}
