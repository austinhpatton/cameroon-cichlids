#!/bin/bash
#SBATCH --output=gather-joint-vcfs_%a.out
#SBATCH --error=gather-joint-vcfs_%a.err
#SBATCH --mail-type=ALL
#SBATCH --time=23:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=austinhpatton@berkeley.edu
#SBATCH --account=ac_fishes
#SBATCH -p savio3

module load java

gatk="/global/home/users/austinhpatton/software/gatk-4.1.8.1/gatk"
BASE=/global/scratch/users/austinhpatton/cichlids/cameroon/Onil_UMD/GATK/PreBQSR/PerScaff/Joint

####
# Bring everything together after joint genotyping. 

/global/home/users/austinhpatton/software/gatk-4.1.8.1/gatk --java-options "-Xmx90G" GatherVcfs \
        -I ${BASE}/NC_031967.2-joint-gt-BM-cichlids-post-bqsr.vcf.gz \
        -I ${BASE}/NC_031972.2-joint-gt-BM-cichlids-post-bqsr.vcf.gz \
        -I ${BASE}/NC_031986.2-joint-gt-BM-cichlids-post-bqsr.vcf.gz \
        -I ${BASE}/NC_031971.2-joint-gt-BM-cichlids-post-bqsr.vcf.gz \
        -I ${BASE}/NC_031965.2-joint-gt-BM-cichlids-post-bqsr.vcf.gz \
        -I ${BASE}/NC_031979.2-joint-gt-BM-cichlids-post-bqsr.vcf.gz \
        -I ${BASE}/NC_031970.2-joint-gt-BM-cichlids-post-bqsr.vcf.gz \
        -I ${BASE}/NC_031980.2-joint-gt-BM-cichlids-post-bqsr.vcf.gz \
        -I ${BASE}/NC_031976.2-joint-gt-BM-cichlids-post-bqsr.vcf.gz \
        -I ${BASE}/NC_031985.2-joint-gt-BM-cichlids-post-bqsr.vcf.gz \
        -I ${BASE}/NC_031981.2-joint-gt-BM-cichlids-post-bqsr.vcf.gz \
        -I ${BASE}/NC_031982.2-joint-gt-BM-cichlids-post-bqsr.vcf.gz \
        -I ${BASE}/NC_031977.2-joint-gt-BM-cichlids-post-bqsr.vcf.gz \
        -I ${BASE}/NC_031984.2-joint-gt-BM-cichlids-post-bqsr.vcf.gz \
        -I ${BASE}/NC_031966.2-joint-gt-BM-cichlids-post-bqsr.vcf.gz \
        -I ${BASE}/NC_031987.2-joint-gt-BM-cichlids-post-bqsr.vcf.gz \
        -I ${BASE}/NC_031974.2-joint-gt-BM-cichlids-post-bqsr.vcf.gz \
        -I ${BASE}/NC_031969.2-joint-gt-BM-cichlids-post-bqsr.vcf.gz \
        -I ${BASE}/NC_031978.2-joint-gt-BM-cichlids-post-bqsr.vcf.gz \
        -I ${BASE}/NC_031975.2-joint-gt-BM-cichlids-post-bqsr.vcf.gz \
        -I ${BASE}/NC_031983.2-joint-gt-BM-cichlids-post-bqsr.vcf.gz \
        -I ${BASE}/NC_031973.2-joint-gt-BM-cichlids-post-bqsr.vcf.gz \
        -I ${BASE}/ShortScaffs-joint-gt-BM-cichlids-post-bqsr.vcf.gz \
        --CREATE_INDEX false \
        -O ${BASE}/all-BM-cichlids-joint-gt-post-bqsr-RAW.vcf

~/software/htslib-1.11/bgzip ${BASE}/all-BM-cichlids-joint-gt-post-bqsr-RAW.vcf
~/software/htslib-1.11/tabix -p vcf ${BASE}/all-BM-cichlids-joint-gt-post-bqsr-RAW.vcf.gz
