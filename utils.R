##-----Define function----
myFileName <- function(prefix,suffix){
  res <- paste0(prefix,"_",format(Sys.time(),"%Y%m%d"),suffix)
  return(res)
}

MyHeatmapSave <- function(ph,prefix,width,height,legend.position="right"){
    pdf(file = myFileName(prefix = prefix,suffix = ".pdf"),
        width = width,height = height)
    draw(ph,
         heatmap_legend_side = legend.position)
    dev.off()
    
    jpeg(file = myFileName(prefix = prefix,suffix = ".jpg"),
         width = width,height = height,units = "in",res = 350)
    draw(ph,heatmap_legend_side = legend.position)
    dev.off()
  
}

calc_ht_size = function(ht, unit = "inch") {
  pdf(NULL)
  ht = draw(ht)
  w = ComplexHeatmap:::width(ht)
  w = convertX(w, unit, valueOnly = TRUE)
  h = ComplexHeatmap:::height(ht)
  h = convertY(h, unit, valueOnly = TRUE)
  dev.off()
  c(w, h)
}




###------ define colors ----------------

rdbu <- colorRampPalette(rev(RColorBrewer::brewer.pal(11,"RdBu"))) 
test.color.2 <- colorRampPalette(c("orange","lightgoldenrod","darkgreen"))
col.spectral <- colorRampPalette(brewer.pal(11,'Spectral')[-6])
test.color.3 <- colorRampPalette(c("#f86e11","#e9970a","#71a701","#62b474","#71c3ac","#9fc4ca"))
rdwhbu <- colorRampPalette(c("navy", "white", "brown3"))
skyblueYelow <- colorRampPalette(c("skyblue","black","yellow"))
skybluered <- colorRampPalette(c("skyblue","black","orange"))
solarExtra <- colorRampPalette(c("#3361A5","#248AF3","#14B3FF","#88CEEF","#C1D5DC","#EAD397","#FDB31A","#E42A2A","#A31D1D"))
blues <- colorRampPalette(colors = brewer.pal(9,"Blues"))
ylord <- colorRampPalette(colors = brewer.pal(9,"YlOrRd"))
hic.red <- colorRampPalette(c("white","red"))
# "#FFF5EB" "#FEE6CE" "#FDD0A2" "#FDAE6B" "#FD8D3C" "#F16913" "#D94801" "#A63603"
# "#7F2704"
hic.orage <- colorRampPalette(c("#FFF5EB","#EAD397","#FDB31A","#E42A2A","#A31D1D"))
sc.hic.orange <- colorRampPalette(c("grey85","#EAD397","#FDB31A","#E42A2A","#A31D1D"))

cold <- colorRampPalette(c('#f7fcf0','#41b6c4','#253494','#081d58'))
warm <- colorRampPalette(c('#ffffb2','#fecc5c','#e31a1c','#800026'))
mypalette <- c(rev(cold(21)), warm(20))
coldwarm <- colorRampPalette(colors = mypalette)

onlywarm <-colorRampPalette(colors = c("#f7fcf0", "#E5D8D1", "#F2CBB7",
                                       "#F7B89C", "#F6A081", "#EE8568",
                                       "#E0654F", "#CC4039", "#B40426"))
coolwarm <- colorRampPalette(colors = c("#3B4CC0", "#4F6AD9", "#6585EC", "#7B9FF9", "#93B5FF", 
                                        "#AAC7FD", "#C0D4F5", "#D4DBE6", "#E5D8D1", "#F2CBB7",
                                        "#F7B89C", "#F6A081", "#EE8568",
                                        "#E0654F", "#CC4039", "#B40426"))
ramp <- colorRampPalette(c("white","pink","red","black"))

bkcolor <- c(colorRampPalette(c(brewer.pal(9,"Blues" )[4:9],"#1a1919"))(50),
             colorRampPalette(c("#1a1919",rev(brewer.pal( 9,"YlOrBr" ))[1:6]))(50))

bkcolor <- colorRampPalette(colors = bkcolor)
hic.pca.red <- colorRampPalette(c("blue","gray1","red"))
hic.pca.redwhite <- colorRampPalette(c("#1d1856","navyblue","white","red4","#861617"))
hic.pca.orange <- colorRampPalette(c("#2f2583","black","#f9b232"))
hic.pca.skyblue <- colorRampPalette(c("skyblue","black","orange"))

okabe_ito <- colorRampPalette(colors = c("#e59f01","#56b4e8","#009f73","#f0e442","#0072b1","#d55e00","#cc79a7","#999999","#000000"))
OrBl_div <- colorRampPalette(colors = c("#9f3d22","#be4d21","#db6525","#ef8531","#f1ac73","#d8d4c9","#a1bccf","#6fa3cb","#5689b6","#4171a1","#2b5b8b"))
OrBl_v2 <- colorRampPalette(colors = rev(c("#940025","#e50c2f","#F15F30" ,"#F7962E", "#FCEE2B",
                                       "#88CEEF","#248AF3","#14B3FF","#3361A5","#004377")))

light_molande <- colorRampPalette(colors = c("#413C58","#A3C4BC","#BFD7B5","#E7EFC5","#F2DDA4"))


coldwarm2.0 <- colorRampPalette(c("#3C1518","#69140E","#A44200","#D58936","#F2F3AE",
                                  "#D4E4BC","#96ACB7","#36558F","#40376E","#48233C"))

###gundam
strike_freedom_gundam_color <- colorRampPalette(colors = c("brown3","#eff3ff","navy","#373b35","#b1b0c2","#f2be58"))
npg.color <- colorRampPalette(colors = ggsci:::ggsci_db$npg$nrc) 

my_stepped3.color <- colorRampPalette(colors = pals::stepped3(20))
my_stepped.color <- colorRampPalette(colors = pals::stepped(24))
#mycolor.bar(hic.pca.red(100),min = -1,max = 1)
divergentcolor <- function (n) {
  colorSpace <- c("#E41A1C", "#377EB8", "#4DAF4A", 
                  "#984EA3", "#F29403", "#F781BF", "#BC9DCC", 
                  "#A65628", "#54B0E4", "#222F75", "#1B9E77", 
                  "#B2DF8A", "#E3BE00", "#FB9A99", "#E7298A", 
                  "#910241", "#00CDD1", "#A6CEE3", "#CE1261", 
                  "#5E4FA2", "#8CA77B", "#00441B", "#DEDC00", 
                  "#B3DE69", "#8DD3C7", "#999999")
  if (n <= length(colorSpace)) {
    colors <- colorSpace[1:n]
  }
  else {
    colors <- (grDevices::colorRampPalette(colorSpace))(n)
  }
  return(colors)
}

####some more divergent color
####from monet
monet.sunset <- colorRampPalette(colors = c("#272924","#323c48","#67879c",
                                            "#7d9390","#d7695a","#ba9f84"))
monet.sunumbrella.women <- colorRampPalette(colors = c("#aed1e4","#bcdde2","#e8e7d5",
                                                       "#f1ede9","#e6cec2","#e49d40",
                                                       "#ddcb3e","#a39524","#89af7a",
                                                       "#744054","#aea3a9"))
monet.cliff <- colorRampPalette(colors = c("#242a26","#3c5653","#bacbc3","#3d4b6e",
                                           "#6f799d","#b3bedc","#f0f0f2","#898788",
                                           "#ddd0bd"))
####show pictures
pic.green <- colorRampPalette(colors = c("#344A23","#726342","#383812","#879979"))
pic.blue <- colorRampPalette(colors = c("#0B2C26","#004054","#95B7B6","#DFE5E1","#C26441"))
pic.green2 <- colorRampPalette(colors = c("#61776B","#A6AA8B","#B6C7D7","#EFEBEC","#E79A61"))
pic.beauty <- colorRampPalette(colors = c("#2F5365","#8A8E8F","#435D2E","#09200B","#CDBBA4","#D04A49"))
pic.beauty2 <- colorRampPalette(colors = c("#4B6032","#436A4D","#DFB28D","#BF4E31","#87988E"))
pic.beauty3 <- colorRampPalette(colors = c("#DB9371","#FCE6D8","#EB9B63","#DBD6C9","#666E20"))
pic.greenblue <- colorRampPalette(colors = c("#86AFA9","#E9A419","#FC8550","#F7C9BA","#8C5E51"))
pic.boy <- colorRampPalette(colors = c("#7C7F62","#956A3D","#625A58","#B1B1B1","#D8D4D5"))
pic.fellows <- colorRampPalette(colors = c("#193B14","#85A063","#8A4925","#C5947E","#ACA9A2"))
pic.girl <- colorRampPalette(colors = c("#1B3E64","#BE1705","#D77363","#84865D","#E2B88D"))
pic.girl2 <- colorRampPalette(colors = c("#BF593F","#DAAE7D","#AAC5DA","#DAC1BA","#141A22"))
molandi.color <- c("#e4d4c5","#c6b1ac","#764e56",
                   "#f9d9c0","#d19477","#93675a",
                   "#f0e9cd","#b9a783","#796656",
                   "#cdc1b1","#a2967e","#656356",
                   "#d8e7e4","#9eb2b1","#5a6873",
                   "#ccd8b0","#7e8563","#50463d")
molandi.color <- colorRampPalette(colors = molandi.color)
#scales::show_col(pic.girl2(10))
pinkblack<- colorRampPalette(colors = c("#333333","#666A86","#95B8D1","#E8DDB5","#EDAFB8"))

beaches <- colorRampPalette(colors = c("#87D2DB" ,"#5BB1CB", "#4F66AF" ,"#F15F30" ,"#F7962E", "#FCEE2B" ))

molandi.color3 <- colorRampPalette(colors = c("#FCD0A1","#B1B695","#A690A4","#5E4B56","#AFD2E9"))

sun_gradual <- colorRampPalette(colors = c("#39489f", "#39bbec", "#f9ed36", "#f38466", "#b81f25"))


rdbu_dp <- colorRampPalette(colors = c("#004377","#007ebe","#00a7d1","#81d1e6","#d0ebf4",
                                       "#f9f9f9","#ffe1cf","#ffb18f","#ff6e58","#e50c2f","#940025"))

# tmp.color <- colorRampPalette(colors = c("#612C6D","#D46889","#B2B2B2",
#                                          "#869CC3","#895C48"))
# 
# tmp.color <- c("#73808E","#86BDD8","#fae79b","#F46F43",
#                "#edb386","#344F99","#895C48","#DCD7C1")

# tmp.color <- c("#51c4c2","#0d8a8c","#4583b3","#f78e26","#f172ad","#f7afb9",
#                "#c63596","#be86ba","#8b66b8","#4068b2","#512a93","#223271")

science.color <- colorRampPalette(c("#b03d26","#005f81",
                                    "#9ccfe6","#e0897e",
                                    "#a5a7ab","black"))



YlGnBu_gradient <- colorRampPalette(rev(hcl.colors(100,palette = "YlGnBu")))

YlGnBu_gradientv2 <- colorRampPalette(colors = brewer.pal(n = 9,name = "YlGnBu"))


BlYl_gradient <- colorRampPalette(colors = rev(hcl.colors(100,palette = "Blue-Yellow")))

#tmp.color <- c('#000000', '#380000', '#560000', '#760100', '#980300', '#bb0600', '#df0d00', '#f93500', '#fe6800', '#ff9100', '#ffb402', '#ffd407', '#fff324')
tmp.color <- c("#000000","#110000","#230000","#350000",'#980300',"#EB4F00","#FFC705","#FFEC1E","#FFEE9D")

myfirev2.color <- colorRampPalette(tmp.color)

tmp.color <- c('#000000', '#380000', '#560000', '#760100', '#980300', '#bb0600', '#df0d00', '#f93500', '#fe6800', '#ff9100', '#ffb402', '#ffd407','#fff324','#ffee9d','#fdf7d2')
myfire.color <- colorRampPalette(tmp.color[c(1,2,7,11,12,13,14)])
