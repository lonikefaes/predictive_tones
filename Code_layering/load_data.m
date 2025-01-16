function [D] = load_data(file)

Info   = spm_vol(file);
D.data = spm_read_vols(Info);
D.info = Info;
D.file = file;