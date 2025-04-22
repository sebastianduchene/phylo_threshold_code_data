import os, re

os.chdir("/home/sebastiand/Dropbox/projects_WORKING/temporal_signal_comment/phylo_threshold_code_data/sampling_window_width/reasonable_prior")

os.chdir("/home/sebastiand/Dropbox/projects_WORKING/temporal_signal_comment/phylo_threshold_code_data/sampling_window_width/misleading_prior")


fasta_files = [f for f in os.listdir(".") if len(re.findall("fasta", f)) > 0]
fasta_files

for f in fasta_files:
    data_set_name = re.sub(".fasta", "", f)
    print(data_set_name)
    if "under_threshold" in data_set_name:
        pass
        # Change to at_threshold for the fasta, xml, and sample_prior.xml. Also change name in xml files
        new_data_set_name = re.sub("_under_threshold", "_at_threshold", data_set_name)
        sample_posterior = open(data_set_name+"_misleading_priors.xml", "r").read()
        sample_posterior = re.sub("_under_threshold", "_at_threshold", sample_posterior)
        open(new_data_set_name+"_misleading_priors.xml", "w").write(sample_posterior)

        sample_prior = open(data_set_name+"_sample_prior_misleading_priors.xml", "r").read()
        sample_prior = re.sub("_under_threshold", "_at_threshold", sample_prior)
        open(new_data_set_name+"_sample_prior_misleading_priors.xml", "w").write(sample_prior)

        alignment = open(data_set_name+".fasta", "r").read()
        open(new_data_set_name+".fasta", "w").write(alignment)
        print(f)
    elif "0.5X_pd_threshold" in data_set_name:
        # this is OK
#        print(f)
        pass
    elif("2X_pd_threshold" in data_set_name):
        # This is actually 10X
        new_data_set_name = re.sub("_2X_", "_10X_", data_set_name)
        sample_posterior = open(data_set_name+"_misleading_priors.xml", "r").read()
        sample_posterior = re.sub("_2X_", "_10X_", sample_posterior)
        open(new_data_set_name+"_misleading_priors.xml", "w").write(sample_posterior)

        sample_prior = open(data_set_name+"_sample_prior_misleading_priors.xml", "r").read()
        sample_prior = re.sub("_2X_", "_10X_", sample_prior)
        open(new_data_set_name+"_sample_prior_misleading_priors.xml", "w").write(sample_prior)

        alignment = open(data_set_name+".fasta", "r").read()
        open(new_data_set_name+".fasta", "w").write(alignment)
        print(f)
    elif("10X_pd_threshold" in data_set_name):
        # This is actually 100X
#        new_data_set_name = re.sub("_10X_", "_100X_", data_set_name)
#        sample_posterior = open(data_set_name+"_misleading_priors.xml", "r").read()
#        sample_posterior = re.sub("_10X_", "_100X_", sample_posterior)
#        open(new_data_set_name+"_misleading_priors.xml", "w").write(sample_posterior)

#        sample_prior = open(data_set_name+"_sample_prior_misleading_priors.xml", "r").read()
#        sample_prior = re.sub("_10X_", "_100X_", sample_prior)
#        open(new_data_set_name+"_sample_prior_misleading_priors.xml", "w").write(sample_prior)

#        alignment = open(data_set_name+".fasta", "r").read()
#        open(new_data_set_name+".fasta", "w").write(alignment)
        #print(f)
        pass
