# ----
library(msa)
library(Biostrings)


ICAM1 <- readAAStringSet("ICAM1.txt")
ICAM1


Organisms <- c("Homo sapiens", "Papio anubis", "Nomascus leucogenys", "Saimiri boliviensis boliviensis",
               "Aotus nancymaae", "Callithrix jacchus", "Neomonachus schauinslandi", "Microcebus murinus",
               "Balaenoptera musculus", "Otolemur garnettii", "Propithecus coquereli", "Ailuropoda melanoleuca",
               "Tupaia chinensis", "Pteropus vampyrus", "Lynx pardinus", "Myotis brandtii", "Lagenorhynchus obliquidens",
               "Marmota marmota", "Galeopterus variegatus", "Hipposideros armiger",
               "Mustela putorius furo", "Choloepus didactylus", "Manis pentadactyla", "Suricata suricatta", "Neophocaena asiaeorientalis",
               "Callorhinus ursinus", "Oryctolagus cuniculus", "Lipotes vexillifer", "Physeter catodon", "Canis lupus familiaris",
               "Bos taurus", "Sus scrofa", "Equus caballus", "Loxodonta africana", "p53")


length(Organisms)


names(ICAM1) <- Organisms

ICAM1_2 <-msa(ICAM1, method = "ClustalW")

library(seqinr)

ICAM1_3 <- msaConvert(ICAM1_2, type="seqinr::alignment")

d <- dist.alignment(ICAM1_3, "identity")
as.matrix(d)

library(ape)

ICAM1_tree <- nj(d)

tre <- root(ICAM1_tree, 35, resolve.root = T)
is.rooted(tre)

plot(tre, main="ICAM1")





# install.packages("TreeTools")

library(TreeTools)

# as.Newick(tre, file = "Icam1_arbol") esto no funciona

# write.tree(phy = tre, file="ICAM1_tree.newick")
 
#----





# Red D1

D1_aln <- readAAMultipleAlignment(file = "D1.1.fa", format = "fasta")
D1_aln

rownames(D1_aln) <- c("Balaenoptera musculus", "Physeter catodon", "Lagenorhynchus obliquidens", "Neophocaena asiaeorientalis", "Lipotes vexillifer",
                      "Bos taurus", "Neomonachus schauinslandi", "Callorhinus ursinus", "Ailuropoda melanoleuca",
                      "Mustela putorius", "Canis lupus familiaris", "Lynx pardinus", "Suricata suricatta", "Manis pentadactyla",
                      "Homo sapiens", "Nomascus leucogenys", "Papio anubis", "Saimiri boliviensis", "Callithrix jacchus",
                      "Aotus nancymaae", "Microcebus murinus",  "Propithecus coquereli", "Otolemur garnettii", "Tupaia chinensis",
                      "Marmota marmota", "Galeopterus variegatus", "Oryctolagus cuniculus", "Choloepus didactylus", "Pteropus vampyrus",
                      "Hipposideros armiger", "Myotis brandtii", "Loxodonta africana", "Equus caballus", "Sus scrofa")

length(rownames(D1_aln))


D1 <- msaConvert(D1_aln, type="seqinr::alignment")

D1 <- dist.alignment(D1, "identity")
D1 <- as.matrix(D1)

class(D1)

View(D1)

heatmap(D1)

library("lattice")
levelplot(D1)


library(pheatmap)

pheatmap(D1, main = "ICAM-1 Dominio 1", clustering_method = "complete", cluster_rows = FALSE, 
         cluster_cols = FALSE, color = colorRampPalette (c ("coral3", "lightyellow", "darkslateblue")) (50))



# Dendograma con hclust
hc <- hclust(as.dist(D1), method = "complete")
plot(hc, cex = 0.6, hang = -1, main = "ICAM-1 Dominio 1")

tre2 <- as.phylo(hc)

# write.tree(phy = tre2, file="D1_dendo.newick")


?hclust




## Tablas
Orden <- c("Primates", 'Primates', 'Primates', 'Primates', 'Primates',
           'Primates', 'Carnivora', 'Primates', 'Artiodactyla', 'Primates',
           'Primates', 'Carnivora', 'Scandentia', 'Chiroptera', 'Carnivora',
           'Chiroptera', 'Artiodactyla', 'Rodentia', 'Dermoptera', 'Chiroptera',
           'Carnivora', 'Pilosa', 'Pholidota', 'Carnivora', 'Artiodactyla',
           'Carnivora', 'Lagomorpha', 'Artiodactyla', 'Artiodactyla', 'Carnivora',
           'Artiodactyla', 'Artiodactyla', 'Perissodactyla', 'Proboscidea')

Familia <- c("Hominidae", "Cercopithecidae", "Hylobatidae", "Cebidae", "Aotidae", 
            "Callitrichidae", "Phocidae", "Cheirogaleidae", "Balaenopteridae", "Galagidae",
            "Indriidae", "Ursidae", "Tupaiidae", "Pteropodidae", "Felidae",
            "Vespertilionidae", "Delphinidae", "Sciuridae", "Cynocephalidae", "Hipposideridae",
            "Mustelidae", "Choloepodidae", "Manidae", "Herpestidae", "Phocoenidae",
            "Otariidae", "Leporidae", "Iniidae", "Physeteridae", "Canidae",
            "Bovidae", "Suidae", "Equidae", "Elephantidae")

Nombre_cientifico <- c("Homo sapiens", "Papio anubis", "Nomascus leucogenys", "Saimiri boliviensis boliviensis","Aotus nancymaae", 
                "Callithrix jacchus", "Neomonachus schauinslandi", "Microcebus murinus","Balaenoptera musculus", "Otolemur garnettii", 
                "Propithecus coquereli", "Ailuropoda melanoleuca", "Tupaia chinensis", "Pteropus vampyrus", "Lynx pardinus", 
                "Myotis brandtii", "Lagenorhynchus obliquidens", "Marmota marmota", "Galeopterus variegatus", "Hipposideros armiger", 
                "Mustela putorius furo", "Choloepus didactylus", "Manis pentadactyla", "Suricata suricatta", 
                "Neophocaena asiaeorientalis", "Callorhinus ursinus", "Oryctolagus cuniculus", "Lipotes vexillifer", "Physeter catodon", 
                "Canis lupus familiaris", "Bos taurus", "Sus scrofa", "Equus caballus", "Loxodonta africana")

Nombre_comun <- c('Humano', 'Babuino de Anubis', 'Gibón de mejillas blancas del norte', 'Mono ardilla boliviano', 'Mono nocturno de Nancy Ma',
                  'Tití común', 'Foca monje de Hawái', 'Lémur ratón gris', 'Ballena azul', 'Gálago de Garnet',
                  'Sifaca de Coquerel', 'Panda', 'Tupaya', 'Gran zorro volador', 'Lince ibérico',
                  'Murciélago', 'Delfín de costados blancos', 'Marmota alpina', 'Kaguang', 'Murciélago de nariz de hoja',
                  'Hurón', 'Perezoso de dos dedos de Linnaeus', 'Pangolín chino', 'Suricata',
                  'Marsopa lisa', 'Oso marino ártico', 'Conejo', 'Baiji', 'Cachalote',
                  'Perro', 'Vaca/Toro', 'Jabalí', 'Caballo', 'Elefante africano de sabana')



Tab_info <- cbind(Familia, Nombre_cientifico, Nombre_comun, Orden)

# write.csv(Tab_info, file = "Info_organismos.csv")






# Idea que no sirvio
#----
range(d2)

D1_aln2[D1_aln2 < .7] <- 0
D1_aln2[D1_aln2 > .69999] <- 1
D1_aln2

library(igraph)

g <- graph.adjacency(as.matrix(D1_aln2), weighted = TRUE, mode = "undirected")
plot(g)

layout <- layout.fruchterman.reingold(g)
layout

plot(g, layout = layout, vertex.size = 8)

plot(degree_distribution(g), main = "Degree distribution", xlab = "Degree", ylab = "Frequency")


table(d2)
