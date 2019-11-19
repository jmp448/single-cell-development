library(stringr)

for (i in 1:3) {
  for (j in 1:6) {

    demux <- read.table(paste0("/project2/gilad/reem/singlecellCM/round1/fulldata/CD",
      i, "/CD", i, "col", j, "/demux/hpCD", i, "col", j, "_demux.best"), header = T,
      stringsAsFactors = F)

    assign(paste0("CD", i, "col", j, "_demux"), demux)
  }
}

save.image(file="./workspaces/demuxlet_data.RData")
