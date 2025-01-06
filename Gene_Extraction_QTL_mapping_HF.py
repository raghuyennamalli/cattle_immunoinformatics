
# Variation:Substitution, Insertions, Deletions
def indicineprotvar(hfprotvar_wd, chr):
  #print("hffinal_df", hfprotvar_wd.shape)
  in_f1=f"ind_chr{chr}_prot.csv"
  ind_prot=pd.read_csv(in_f1)
  ind_prot["indstart_of_next_gene"] = ind_prot["Start"].shift(-1).fillna(0)
  finalindprot_df=pd.DataFrame(columns=list(hfprotvar_wd)+["NelStart"]+["NelStop"]+["NelLocus"]+["NelLength"]+["NelProt"]+["Nelnextgen"]+["Nelgene_type"])
  for i, row in tqdm(ind_prot.iterrows()):
    # Check for the at start
    if hfprotvar_wd[(row["Start"]==hfprotvar_wd["startpos.q"]) & (row["Stop"]>hfprotvar_wd["endpos.q"])].shape[0]>=1:
        tmp_df = hfprotvar_wd[(row["Start"]==hfprotvar_wd["startpos.q"]) & (row["Stop"]>hfprotvar_wd["endpos.q"])].copy(deep=True)
        tmp_df["Nelgene_type"] = "at_start"
        tmp_df["NelStart"] = row["Start"]
        tmp_df["NelStop"] = row["Stop"]
        tmp_df["NelLocus"] = row["Locus"]
        tmp_df["NelLength"] = row["Length"]
        tmp_df["NelProt"] = row["Protein Name"]
        tmp_df["Nelnextgen"] = row["indstart_of_next_gene"]
        finalindprot_df = pd.concat([finalindprot_df,tmp_df]).reset_index(drop=True)
    # check for the at end
    elif hfprotvar_wd[(row["Start"]<hfprotvar_wd["startpos.q"]) & (row["Stop"]==hfprotvar_wd["endpos.q"])].shape[0]>=1:
        tmp_df = hfprotvar_wd[(row["Start"]<hfprotvar_wd["startpos.q"]) & (row["Stop"]==hfprotvar_wd["endpos.q"])].copy(deep=True)
        tmp_df["Nelgene_type"] = "at_end"
        tmp_df["NelStart"] = row["Start"]
        tmp_df["NelStop"] = row["Stop"]
        tmp_df["NelLocus"] = row["Locus"]
        tmp_df["NelLength"] = row["Length"]
        tmp_df["NelProt"] = row["Protein Name"]
        tmp_df["Nelnextgen"] = row["indstart_of_next_gene"]
        finalindprot_df = pd.concat([finalindprot_df,tmp_df]).reset_index(drop=True)
    # Check for the within gene
    elif hfprotvar_wd[(row["Start"]<hfprotvar_wd["startpos.q"]) & (row["Stop"]>hfprotvar_wd["endpos.q"])].shape[0]>=1:
        tmp_df = hfprotvar_wd[(row["Start"]<hfprotvar_wd["startpos.q"]) & (row["Stop"]>hfprotvar_wd["endpos.q"])].copy(deep=True)
        tmp_df["Nelgene_type"] = "within_gene"
        tmp_df["NelStart"] = row["Start"]
        tmp_df["NelStop"] = row["Stop"]
        tmp_df["NelLocus"] = row["Locus"]
        tmp_df["NelLength"] = row["Length"]
        tmp_df["NelProt"] = row["Protein Name"]
        tmp_df["Nelnextgen"] = row["indstart_of_next_gene"]
        finalindprot_df = pd.concat([finalindprot_df,tmp_df]).reset_index(drop=True)
    # Check for the intergenic gene
    elif hfprotvar_wd[(row["Stop"]<hfprotvar_wd["startpos.q"]) & (row["indstart_of_next_gene"]>hfprotvar_wd["endpos.q"])].shape[0]>=1:
        tmp_df = hfprotvar_wd[(row["Stop"]<hfprotvar_wd["startpos.q"]) & (row["indstart_of_next_gene"]>hfprotvar_wd["endpos.q"])].copy(deep=True)
        tmp_df["Nelgene_type"] = "intergenic_gene"
        tmp_df["NelStart"] = row["Start"]
        tmp_df["NelStop"] = row["Stop"]
        tmp_df["NelLocus"] = row["Locus"]
        tmp_df["NelLength"] = row["Length"]
        tmp_df["NelProt"] = row["Protein Name"]
        tmp_df["Nelnextgen"] = row["indstart_of_next_gene"]
        finalindprot_df = pd.concat([finalindprot_df,tmp_df]).reset_index(drop=True)
    else:
        pass


  # Generate output file
  out_file20=f"hfqtl_NelVar{chr}.csv"
  finalindprot_df.to_csv(out_file20)
  #shutil.copy(out_file13, outloc)

  # extract only the immune genes
  hfimm_nel = pd.DataFrame
  hfimm_nel=finalindprot_df[finalindprot_df["keyword_found"]==True]
  out_file21=f"hfimm_nel{chr}.csv"
  hfimm_nel.to_csv(out_file21)
  #shutil.copy(out_file14, outloc)

  # LOV50
  hfNel_LOV50=finalindprot_df[finalindprot_df["LOV"]>=50]
  hfNelimm_LOV50=hfimm_nel[hfimm_nel["LOV"]>=50]
  out_file22=f"hfNel_LOV50_{chr}.csv"
  hfNel_LOV50.to_csv(out_file22)
  out_file23=f"hfNelimm_LOV50_{chr}.csv"
  hfNelimm_LOV50.to_csv(out_file23)
  #

def hfqtl_VarAnal(common_df, hfqtl_df, chr):
  # Identifying the variations
  final_df = pd.DataFrame(columns=list(hfqtl_df.columns)+["HFgene_type"]+["startpos.r"]+["endpos.r"]+["startpos.q"]+["endpos.q"]+["seq.r"]+["seq.q"]+["LOV"])
  for i, row in tqdm(common_df.iterrows()):
    # Check for the at start
    if hfqtl_df[(hfqtl_df["Begin"]==row["startpos.r"]) & (hfqtl_df["End"]>row["endpos.r"])].shape[0]>=1:
        tmp_df = hfqtl_df[(hfqtl_df["Begin"]==row["startpos.r"]) & (hfqtl_df["End"]>row["endpos.r"])].copy(deep=True)
        tmp_df["HFgene_type"] = "at_start"
        tmp_df["startpos.r"] = row["startpos.r"]
        tmp_df["endpos.r"] = row["endpos.r"]
        tmp_df["startpos.q"] = row["startpos.q"]
        tmp_df["endpos.q"] = row["endpos.q"]
        tmp_df["seq.r"] = row["seq.r"]
        tmp_df["seq.q"] = row["seq.q"]
        tmp_df["LOV"] = row["LOV"]
        final_df = pd.concat([final_df,tmp_df]).reset_index(drop=True)
    # check for the at end
    elif hfqtl_df[(hfqtl_df["Begin"]<row["startpos.r"]) & (hfqtl_df["End"]==row["endpos.r"])].shape[0]>=1:
        tmp_df = hfqtl_df[(hfqtl_df["Begin"]<row["startpos.r"]) & (hfqtl_df["End"]==row["endpos.r"])].copy(deep=True)
        tmp_df["HFgene_type"] = "at_end"
        tmp_df["startpos.r"] = row["startpos.r"]
        tmp_df["endpos.r"] = row["endpos.r"]
        tmp_df["startpos.q"] = row["startpos.q"]
        tmp_df["endpos.q"] = row["endpos.q"]
        tmp_df["seq.r"] = row["seq.r"]
        tmp_df["seq.q"] = row["seq.q"]
        tmp_df["LOV"] = row["LOV"]
        final_df = pd.concat([final_df,tmp_df]).reset_index(drop=True)
    # Check for the within gene
    elif hfqtl_df[(hfqtl_df["Begin"]<row["startpos.r"]) & (hfqtl_df["End"]>row["endpos.r"])].shape[0]>=1:
        tmp_df = hfqtl_df[(hfqtl_df["Begin"]<row["startpos.r"]) & (hfqtl_df["End"]>row["endpos.r"])].copy(deep=True)
        tmp_df["HFgene_type"] = "within_gene"
        tmp_df["startpos.r"] = row["startpos.r"]
        tmp_df["endpos.r"] = row["endpos.r"]
        tmp_df["startpos.q"] = row["startpos.q"]
        tmp_df["endpos.q"] = row["endpos.q"]
        tmp_df["seq.r"] = row["seq.r"]
        tmp_df["seq.q"] = row["seq.q"]
        tmp_df["LOV"] = row["LOV"]
        final_df = pd.concat([final_df,tmp_df]).reset_index(drop=True)
    # Check for the intergenic gene
    elif hfqtl_df[(hfqtl_df["End"]<row["startpos.r"]) & (hfqtl_df["hfstart_of_next_gene"]>row["endpos.r"])].shape[0]>=1:
        tmp_df = hfqtl_df[(hfqtl_df["End"]<row["startpos.r"]) & (hfqtl_df["hfstart_of_next_gene"]>row["endpos.r"])].copy(deep=True)
        tmp_df["HFgene_type"] = "intergenic_gene"
        tmp_df["startpos.r"] = row["startpos.r"]
        tmp_df["endpos.r"] = row["endpos.r"]
        tmp_df["startpos.q"] = row["startpos.q"]
        tmp_df["endpos.q"] = row["endpos.q"]
        tmp_df["seq.r"] = row["seq.r"]
        tmp_df["seq.q"] = row["seq.q"]
        tmp_df["LOV"] = row["LOV"]
        final_df = pd.concat([final_df,tmp_df]).reset_index(drop=True)
    else:
        pass

  # Generate output file
  out_file3=f"hfprotvar{chr}.csv"
  final_df.to_csv(out_file3)
  #shutil.copy(out_file1, outloc)

  # find variation not present in gene and intergene region called as other non coding region
  final_df1 = pd.DataFrame(columns=common_df.columns)
  for index, row in tqdm(common_df.iterrows()):
    # Check if 'startpos.q' is present in any row of file2
    if row['startpos.q'] not in final_df['startpos.q'].values:
        # Append the row to the DataFrame if not found
        final_df1 = final_df1.append(row, ignore_index=True)

  # Save the DataFrame with rows not found in file2 to a new CSV file
  out_file4=f"non_coding{chr}.csv"
  final_df1.to_csv(out_file4, index=False)
  #shutil.copy(out_file2, outloc)

  # remove duplicates
  hfprotvar_wd = final_df.drop_duplicates(varet=['Symbol', "QTL Class", 'startpos.r', 'endpos.r'])
  out_file5=f"hfprotvar_wd{chr}.csv"
  hfprotvar_wd.to_csv(out_file5)
  #shutil.copy(out_file_1, outloc)

  # extract only the immune genes
  hfimm = pd.DataFrame
  hfimm=hfprotvar_wd[hfprotvar_wd["keyword_found"]==True]
  out_file6=f"hfimm{chr}.csv"
  hfimm.to_csv(out_file6)
  #shutil.copy(out_file3, outloc)

  # Extract the variation presents in gene and intergene seperately
  hf_within_gene=hfprotvar_wd[hfprotvar_wd["HFgene_type"]=="within_gene"]
  hf_inter_gene=hfprotvar_wd[hfprotvar_wd["HFgene_type"]=="intergenic_gene"]

  # extract the variations among genomic regions including gene and intergene has greater the 50 LOV (length of varaition)
  hf_LOV50=hfprotvar_wd[hfprotvar_wd["LOV"]>=50]
  
  out_file7=f"hf_within_gene{chr}.csv"
  out_file8=f"hf_inter_gene{chr}.csv"
  out_file9=f"hf_LOV50_{chr}.csv"

  hf_within_gene.to_csv(out_file7)
  hf_inter_gene.to_csv(out_file8)
  hf_LOV50.to_csv(out_file9)


  # Extract the variation presents in immune gene and intergene seperately
  hfimm_within_gene=hfimm[hfimm["HFgene_type"]=="within_gene"]
  hfimm_inter_gene=hfimm[hfimm["HFgene_type"]=="intergenic_gene"]

  # extract the variations among genomic regions including gene and intergene has greater the 50 LOV (length of varaition)
  hfimm_LOV50=hfimm[hfimm["LOV"]>=50]

  out_file10=f"hfimm_within_gene{chr}.csv"
  out_file11=f"hfimm_inter_gene{chr}.csv"
  out_file12=f"hfimm_LOV50_{chr}.csv"

  hfimm_within_gene.to_csv(out_file10)
  hfimm_inter_gene.to_csv(out_file11)
  hfimm_LOV50.to_csv(out_file12)

  

  # count the variation based on gene type
  hfvar_counts = hfprotvar_wd.drop_duplicates(varet=['Symbol', 'startpos.r', 'endpos.r'])
  gene_count=hfvar_counts.groupby(["HFgene_type"]).size().reset_index(name="count")
  out_file13=f"gene_count{chr}.csv"
  gene_count.to_csv(out_file13)
  #shutil.copy(out_file10, outloc)


  # extract the immune has how many variations and its gene type
  hfimmvar_counts = hfimm.drop_duplicates(varet=['Symbol', 'startpos.r', 'endpos.r'])
  imm_count=hfimmvar_counts.groupby(by=["Symbol","Name", "HFgene_type"]).size().reset_index(name="count")
  out_file14=f"immlocus_gene_type_counts{chr}.csv"
  imm_count.to_csv(out_file14, index=False)
  #shutil.copy(out_file11, outloc)

  # Calculate value counts
  count = hfprotvar_wd["HFgene_type"].value_counts().reset_index()
  count.columns = ["HFgene_type", "Count"]
  count.loc[len(count.index)]=["immune_gene", len(hfimm_within_gene)]
  count.loc[len(count.index)]=["non_coding_gene", len(final_df1)]


  # Save the result to a CSV file
  out_file15=f"var_count{chr}.csv"
  count.to_csv(out_file15, index=False)
  #shutil.copy(out_file12, outloc)
  
  # count vartitution per gene
  hf_within_gene_counts=hf_within_gene.drop_duplicates(varet=['Symbol', 'startpos.r', 'endpos.r'])
  counts=hf_within_gene_counts.groupby(["Symbol"]).size().reset_index(name="count")
  out_file16=f"loc_var_counts{chr}.csv"
  counts.to_csv(out_file16)
  #shutil.copy(out_file_10, outloc)
  loc50_var=counts[counts["count"]>=50]
  out_file17=f"loc50_var{chr}.csv"
  loc50_var.to_csv(out_file17)
  
  # count vartitution per immune gene
  hfimm_within_gene_counts=hfimm_within_gene.drop_duplicates(varet=['Symbol', 'startpos.r', 'endpos.r'])
  counts=hfimm_within_gene_counts.groupby(["Symbol"]).size().reset_index(name="count")
  out_file18=f"immloc_var_counts{chr}.csv"
  counts.to_csv(out_file18)
  #shutil.copy(out_file_10, outloc)
  immloc50_var=counts[counts["count"]>=50]
  out_file19=f"immloc50_var{chr}.csv"
  immloc50_var.to_csv(out_file19)

  #extract idicine genes
  indicineprotvar(hfprotvar_wd, chr)

def map_keywords(str_sentence):
    list_keywords = [
        "immunoglobulin","immunoreceptor","autoimmune","TLR","IgG",
        "autoimmune","autophagy","immunogen","immune","innate","T-cell","NF-kappa", "antigen",
        "B-cell","lymphocyte","histocompatibility","CD24","CD4","LY96", "BCR",
        "IFIT3","PGLYRP1","NKG2D","UL16","leukocyte","cytokine", "interleukin","interferon"
        "antimicrobial peptide","beta-defensin 2","IL16","IL2","chemokine", "antibody"
    ]
    bool_found = bool(
        [
            i for i in list_keywords if i.lower() in str_sentence.lower()
        ]
    )
    return bool_found




import pandas as pd
from tqdm import tqdm
import shutil

Chr_list=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, "X"]
for chr in Chr_list:
  in_file1=f"common_wd{chr}.csv"
  in_file2=f"HF_QTL{chr}.csv"
  common_df=pd.read_csv(in_file1)
  hfqtl_df=pd.read_csv(in_file2)

  print("common_df", common_df.shape)
  print("hfqtl_df", hfqtl_df.shape)


  common_df["LOV"]=abs(common_df["endpos.q"]-common_df["startpos.q"])
  out_f1=f"common{chr}_LOV.csv"

  common_df.to_csv(out_f1)
  #shutil.copy(outf1, outloc)

  hfqtl_df["keyword_found"] = hfqtl_df.apply(lambda row: map_keywords(row["Name"]),axis=1)
  out_f2=f"hfqtlprot{chr}.csv"

  hfqtl_df.to_csv(out_f2)
  #shutil.copy(outf2, outloc)

  hfqtl_VarAnal(common_df, hfqtl_df, chr)
