library(dplyr)
library(ggplot2)
library(optparse)


# Step 1. Parse arguments
options.list <- list(
    make_option(c("-i", "--tsv-file"),
                type="character",
                dest="tsv.file",
                help="VSTOL TSV file"),
    make_option(c("-o", "--output-pdf-file"),
                type="character",
                dest="output.pdf.file",
                help="Output PDF file")
)
parser <- OptionParser(usage="%prog [options]",
                       option_list = options.list)
args <- parse_args(parser, args = commandArgs(trailingOnly=TRUE))

IsGzipped <- function(file.path) {
  con <- file(file.path, "rb")
  bytes <- readBin(con, "raw", n = 2)
  close(con)
  return(identical(bytes, as.raw(c(0x1f, 0x8b))))
}

RoundUpToNextZeros <- function(x) {
  if (x <= 0) {
    stop("Input must be a positive integer")
  }
  if (x < 10) {
    return(10)
  }
  if (x < 100) {
    return(100)
  }
  magnitude <- 10^(nchar(as.character(x)) - 1)
  rounded <- ceiling(x / magnitude) * magnitude
  if (rounded == x) {
    rounded <- rounded * 10
  }
  return(rounded)
}

# Step 2. Read the variant calls
if (IsGzipped(file.path = args$tsv.file)) {
  df.variant.calls <- read.csv(gzfile(args$tsv.file), sep = "\t", check.names = F, stringsAsFactors = F)
} else {
  df.variant.calls <- read.csv(args$tsv.file, sep = "\t", check.names = F, stringsAsFactors = F)
}

# Plot parameters
PLOT.TITLE.SIZE <- 18
AXIS.TITLE.SIZE <- 16
AXIS.TEXT.SIZE <- 16
BARPLOT.COLOR <- "#000000"
BARPLOT.FILL <- "#D9DEE7"
BARPLOT.COUNT.SIZE <- 4
HISTOGRAM.COLOR <- "#000000"
HISTOGRAM.FILL <- "#D9DEE7"
HISTOGRAM.COUNT.SIZE <- 4

pdf(file = args$output.pdf.file, onefile = TRUE, width = 12, height = 6.75)

# Step 3. Plot the number of variant calls by variant type
df.variant.types <- df.variant.calls %>%
  dplyr::group_by(variant_type) %>%
  summarise(n = n())
plot.variant.types <- ggplot(df.variant.types, aes(x = variant_type, y = n)) +
  geom_bar(stat = "identity", color = BARPLOT.COLOR, fill = BARPLOT.FILL) + 
  xlab("Variant Type") + ylab("Number of Variant Calls") +
  geom_text(stat = 'identity', aes(label = ..y..), vjust = -0.5, size = BARPLOT.COUNT.SIZE) + 
  ggtitle("Number of Variant Calls by variant_type") +
  scale_y_continuous(expand = c(0,1), limits = c(0, RoundUpToNextZeros(max(df.variant.types$n)))) +
  theme_bw() +
  theme(plot.title = element_text(size = PLOT.TITLE.SIZE),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.title = element_text(size = AXIS.TITLE.SIZE),
        axis.text = element_text(size = AXIS.TEXT.SIZE))
print(plot.variant.types)

# Step 4. Plot the number of variant calls by chromosome 1
df.chromosome.1 <- df.variant.calls %>%
  dplyr::group_by(chromosome_1) %>%
  summarise(n = n())
plot.chromosome.1 <- ggplot(df.chromosome.1, aes(x = chromosome_1, y = n)) +
  geom_bar(stat = "identity", color = BARPLOT.COLOR, fill = BARPLOT.FILL) + 
  xlab("Chromosome 1") + ylab("Number of Variant Calls") + 
  ggtitle("Number of Variant Calls by chromosome_1") +
  scale_y_continuous(expand = c(0,1), limits = c(0, RoundUpToNextZeros(max(df.chromosome.1$n)))) +
  geom_text(stat = 'identity', aes(label = ..y..), vjust = -0.5, size = BARPLOT.COUNT.SIZE) + 
  theme_bw() +
  theme(plot.title = element_text(size = PLOT.TITLE.SIZE),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.title = element_text(size = AXIS.TITLE.SIZE),
        axis.text.x = element_text(size = AXIS.TEXT.SIZE, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = AXIS.TEXT.SIZE))
print(plot.chromosome.1)

# Step 5. Plot the number of variant calls by chromosome 2
df.chromosome.2 <- df.variant.calls %>%
  dplyr::group_by(chromosome_2) %>%
  summarise(n = n())
plot.chromosome.2 <- ggplot(df.chromosome.2, aes(x = chromosome_2, y = n)) +
  geom_bar(stat = "identity", color = BARPLOT.COLOR, fill = BARPLOT.FILL) + 
  xlab("Chromosome 2") + ylab("Number of Variant Calls") +
  ggtitle("Number of Variant Calls by chromosome_2") +
  scale_y_continuous(expand = c(0,1), limits = c(0, RoundUpToNextZeros(max(df.chromosome.2$n)))) +
  geom_text(stat = 'identity', aes(label = ..y..), vjust = -0.5, size = BARPLOT.COUNT.SIZE) +
  theme_bw() +
  theme(plot.title = element_text(size = PLOT.TITLE.SIZE),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.title = element_text(size = AXIS.TITLE.SIZE),
        axis.text.x = element_text(size = AXIS.TEXT.SIZE, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = AXIS.TEXT.SIZE))
print(plot.chromosome.2)

# Step 6. Plot the number of variant calls by sample ID
df.sample.id <- df.variant.calls %>%
  dplyr::group_by(sample_id) %>%
  summarise(n = n())
plot.sample.id <- ggplot(df.sample.id, aes(x = sample_id, y = n)) +
  geom_bar(stat = "identity", color = BARPLOT.COLOR, fill = BARPLOT.FILL) + 
  xlab("Sample ID") + ylab("Number of Variant Calls") +
  ggtitle("Number of Variant Calls by sample_id") +
  scale_y_continuous(expand = c(0,1), limits = c(0, RoundUpToNextZeros(max(df.sample.id$n)))) +
  geom_text(stat = 'identity', aes(label = ..y..), vjust = -0.5, size = BARPLOT.COUNT.SIZE) +
  theme_bw() +
  theme(plot.title = element_text(size = PLOT.TITLE.SIZE),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.title = element_text(size = AXIS.TITLE.SIZE),
        axis.text.x = element_text(size = AXIS.TEXT.SIZE, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = AXIS.TEXT.SIZE))
print(plot.sample.id)

# Step 7. Plot the number of variant calls by source ID
df.source.id <- df.variant.calls %>%
  dplyr::group_by(source_id) %>%
  summarise(n = n())
plot.source.id <- ggplot(df.source.id, aes(x = source_id, y = n)) +
  geom_bar(stat = "identity", color = BARPLOT.COLOR, fill = BARPLOT.FILL) + 
  xlab("Source ID") + ylab("Number of Variant Calls") +
  ggtitle("Number of Variant Calls by source_id") +
  scale_y_continuous(expand = c(0,1), limits = c(0, RoundUpToNextZeros(max(df.source.id$n)))) +
  geom_text(stat = 'identity', aes(label = ..y..), vjust = -0.5, size = BARPLOT.COUNT.SIZE) +
  theme_bw() +
  theme(plot.title = element_text(size = PLOT.TITLE.SIZE),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.title = element_text(size = AXIS.TITLE.SIZE),
        axis.text.x = element_text(size = AXIS.TEXT.SIZE, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = AXIS.TEXT.SIZE))
print(plot.source.id)

# Step 8. Plot the number of variant calls by variant calling method
df.variant.calling.method <- df.variant.calls %>%
  dplyr::group_by(variant_calling_method) %>%
  summarise(n = n())
plot.variant.calling.method <- ggplot(df.variant.calling.method, aes(x = variant_calling_method, y = n)) +
  geom_bar(stat = "identity", color = BARPLOT.COLOR, fill = BARPLOT.FILL) + 
  xlab("Variant Calling Method") + ylab("Number of Variant Calls") +
  ggtitle("Number of Variant Calls by variant_calling_method") +
  scale_y_continuous(expand = c(0,1), limits = c(0, RoundUpToNextZeros(max(df.variant.calling.method$n)))) +
  geom_text(stat = 'identity', aes(label = ..y..), vjust = -0.5, size = BARPLOT.COUNT.SIZE) +
  theme_bw() +
  theme(plot.title = element_text(size = PLOT.TITLE.SIZE),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.title = element_text(size = AXIS.TITLE.SIZE),
        axis.text.x = element_text(size = AXIS.TEXT.SIZE, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = AXIS.TEXT.SIZE))
print(plot.variant.calling.method)

# Step 9. Plot the number of variant calls by filter
df.filter <- df.variant.calls %>%
  dplyr::group_by(filter) %>%
  summarise(n = n())
plot.filter <- ggplot(df.filter, aes(x = filter, y = n)) +
  geom_bar(stat = "identity", color = BARPLOT.COLOR, fill = BARPLOT.FILL) + 
  xlab("Filter") + ylab("Number of Variant Calls") +
  ggtitle("Number of Variant Calls by filter") +
  scale_y_continuous(expand = c(0,1), limits = c(0, RoundUpToNextZeros(max(df.filter$n)))) +
  geom_text(stat = 'identity', aes(label = ..y..), vjust = -0.5, size = BARPLOT.COUNT.SIZE) +
  theme_bw() +
  theme(plot.title = element_text(size = PLOT.TITLE.SIZE),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.title = element_text(size = AXIS.TITLE.SIZE),
        axis.text.x = element_text(size = AXIS.TEXT.SIZE, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = AXIS.TEXT.SIZE))
print(plot.filter)

# Step 10. Plot the number of variant calls by precise
df.precise <- df.variant.calls %>%
  dplyr::group_by(precise) %>%
  summarise(n = n())
plot.precise <- ggplot(df.precise, aes(x = precise, y = n)) +
  geom_bar(stat = "identity", color = BARPLOT.COLOR, fill = BARPLOT.FILL) + 
  xlab("Precise") + ylab("Number of Variant Calls") +
  ggtitle("Number of Variant Calls by precise") +
  scale_y_continuous(expand = c(0,1), limits = c(0, RoundUpToNextZeros(max(df.precise$n)))) +
  geom_text(stat = 'identity', aes(label = ..y..), vjust = -0.5, size = BARPLOT.COUNT.SIZE) +
  theme_bw() +
  theme(plot.title = element_text(size = PLOT.TITLE.SIZE),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.title = element_text(size = AXIS.TITLE.SIZE),
        axis.text.x = element_text(size = AXIS.TEXT.SIZE, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = AXIS.TEXT.SIZE))
print(plot.precise)

# Step 11. Plot the number of variant calls by variant size
HelperFn <- function(x) {
  chromosome_1 <- as.character(x[['chromosome_1']])
  chromosome_2 <- as.character(x[['chromosome_2']])
  if (grepl(x[['chromosome_1']], x[['chromosome_2']])) {
    variant.size <- as.integer(x[['variant_size']])
    if (variant.size <= 0) {
      return("Unknown")
    } else if (variant.size >= 1 & variant.size < 50) {
      return("1-50")
    } else if (variant.size >= 50 & variant.size < 100) {
      return("50-100")
    } else if (variant.size >= 100 & variant.size < 500) {
      return("100-500")
    } else if (variant.size >= 500 & variant.size < 1000) {
      return("500-1,000")
    } else if (variant.size >= 1000 & variant.size < 10000) {
      return("1,000-10,000")
    } else if (variant.size >= 10000 & variant.size < 100000) {
      return("10,000-100,000")
    } else if (variant.size >= 100000 & variant.size < 1000000) {
      return("100,000-1,000,000")
    } else {
      return("1,000,000+")
    }
  } else {
    # Interchromosomal
    return("Interchromosomal Event")
  }
}
df.variant.calls.temp <- df.variant.calls
df.variant.calls.temp$variant_size_group <- apply(X = df.variant.calls.temp, MARGIN = 1, FUN = HelperFn)
df.variant.size <- df.variant.calls.temp %>%
  dplyr::group_by(variant_size_group) %>%
  summarise(n = n())
df.variant.size$variant_size_group <- factor(
  x = df.variant.size$variant_size_group,
  levels = c("Unknown","Interchromosomal Event","1-50","50-100","100-500","500-1,000","1,000-10,000","10,000-100,000","100,000-1,000,000","1,000,000+")
)
plot.variant.size <- ggplot(df.variant.size, aes(x = variant_size_group, y = n)) +
  geom_bar(stat = "identity", color = BARPLOT.COLOR, fill = BARPLOT.FILL) + 
  xlab("Variant Size (bp)") + ylab("Number of Variant Calls") +
  ggtitle("Number of Variant Calls by variant_size") +
  scale_y_continuous(expand = c(0,1), limits = c(0, RoundUpToNextZeros(max(df.variant.size$n)))) +
  geom_text(stat = 'identity', aes(label = ..y..), vjust = -0.5, size = BARPLOT.COUNT.SIZE) +
  theme_bw() +
  theme(plot.title = element_text(size = PLOT.TITLE.SIZE),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.title = element_text(size = AXIS.TITLE.SIZE),
        axis.text.x = element_text(size = AXIS.TEXT.SIZE, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = AXIS.TEXT.SIZE))
print(plot.variant.size)

# Step 12. Plot a histogram of the reference allele read count
df.variant.calls.temp <- df.variant.calls[df.variant.calls$reference_allele_read_count >= 0,]
if (nrow(df.variant.calls.temp) > 0) {
  p <- ggplot(df.variant.calls.temp, aes(x = reference_allele_read_count)) + geom_histogram()
  p.data <- ggplot_build(p)$data[[1]]
  max.frequency <- max(p.data$count)
  plot.ref.allele.read.count.histogram <- ggplot(df.variant.calls.temp, aes(x = reference_allele_read_count)) +
    geom_histogram(color = HISTOGRAM.COLOR, fill = HISTOGRAM.FILL, boundary = 0) + 
    xlab("Reference Allele Read Count") + ylab("Frequency") +
    ggtitle("Histogram of reference_allele_read_count") +
    scale_y_continuous(expand = c(0,1), limits = c(0, RoundUpToNextZeros(max.frequency))) +
    theme_bw() +
    theme(plot.title = element_text(size = PLOT.TITLE.SIZE),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.title = element_text(size = AXIS.TITLE.SIZE),
          axis.text = element_text(size = AXIS.TEXT.SIZE))
  print(plot.ref.allele.read.count.histogram)
}

# Step 13. Plot a histogram of the alternate allele read count
df.variant.calls.temp <- df.variant.calls[df.variant.calls$alternate_allele_read_count >= 0,]
if (nrow(df.variant.calls.temp) > 0) {
  p <- ggplot(df.variant.calls.temp, aes(x = alternate_allele_read_count)) + geom_histogram()
  p.data <- ggplot_build(p)$data[[1]]
  max.frequency <- max(p.data$count)
  plot.alt.allele.read.count.histogram <- ggplot(df.variant.calls.temp, aes(x = alternate_allele_read_count)) +
    geom_histogram(color = HISTOGRAM.COLOR, fill = HISTOGRAM.FILL, boundary = 0) + 
    xlab("Alternate Allele Read Count") + ylab("Frequency") +
    ggtitle("Histogram of alternate_allele_read_count") +
    scale_y_continuous(expand = c(0,1), limits = c(0, RoundUpToNextZeros(max.frequency))) +
    theme_bw() +
    theme(plot.title = element_text(size = PLOT.TITLE.SIZE),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.title = element_text(size = AXIS.TITLE.SIZE),
          axis.text = element_text(size = AXIS.TEXT.SIZE))
  print(plot.alt.allele.read.count.histogram) 
}

# Step 14. Plot a histogram of the total read count
df.variant.calls.temp <- df.variant.calls[df.variant.calls$total_read_count >= 0,]
if (nrow(df.variant.calls.temp) > 0) {
  p <- ggplot(df.variant.calls.temp, aes(x = total_read_count)) + geom_histogram()
  p.data <- ggplot_build(p)$data[[1]]
  max.frequency <- max(p.data$count)
  plot.total.read.count.histogram <- ggplot(df.variant.calls.temp, aes(x = total_read_count)) +
    geom_histogram(color = HISTOGRAM.COLOR, fill = HISTOGRAM.FILL, boundary = 0) + 
    xlab("Total Read Count") + ylab("Frequency") +
    ggtitle("Histogram of total_read_count") +
    scale_y_continuous(expand = c(0,1), limits = c(0, RoundUpToNextZeros(max.frequency))) +
    theme_bw() +
    theme(plot.title = element_text(size = PLOT.TITLE.SIZE),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.title = element_text(size = AXIS.TITLE.SIZE),
          axis.text = element_text(size = AXIS.TEXT.SIZE))
  print(plot.total.read.count.histogram)
}

# Step 15. Plot a histogram of the alternate allele fraction
df.variant.calls.temp <- df.variant.calls[df.variant.calls$alternate_allele_fraction > 0.0,]
if (nrow(df.variant.calls.temp) > 0) {
  p <- ggplot(df.variant.calls.temp, aes(x = alternate_allele_fraction)) + geom_histogram(binwidth = 0.05)
  p.data <- ggplot_build(p)$data[[1]]
  max.frequency <- max(p.data$count)
  plot.vaf.histogram <- ggplot(df.variant.calls.temp, aes(x = alternate_allele_fraction)) +
    geom_histogram(binwidth = 0.05, color = HISTOGRAM.COLOR, fill = HISTOGRAM.FILL, boundary = 0) + 
    xlab("Alternate Allele Fraction") + ylab("Frequency") +
    ggtitle("Histogram of alternate_allele_fraction") +
    coord_cartesian(xlim = c(0, 1)) +
    scale_y_continuous(expand = c(0,1), limits = c(0, RoundUpToNextZeros(max.frequency))) +
    theme_bw() +
    theme(plot.title = element_text(size = PLOT.TITLE.SIZE),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.title = element_text(size = AXIS.TITLE.SIZE),
          axis.text = element_text(size = AXIS.TEXT.SIZE))
  print(plot.vaf.histogram)
} 

# Step 16. Plot a histogram of the quality score
df.variant.calls.temp <- df.variant.calls[
  (df.variant.calls$quality_score >= 0) &
  (is.na(df.variant.calls$quality_score) == FALSE),
]
if (nrow(df.variant.calls.temp) > 0) {
  p <- ggplot(df.variant.calls.temp, aes(x = quality_score)) + geom_histogram()
  p.data <- ggplot_build(p)$data[[1]]
  max.frequency <- max(p.data$count)
  plot.quality.score <- ggplot(df.variant.calls.temp, aes(x = quality_score)) +
    geom_histogram(color = HISTOGRAM.COLOR, fill = HISTOGRAM.FILL, boundary = 0) + 
    xlab("Quality Score") + ylab("Frequency") +
    ggtitle("Histogram of quality_score") +
    scale_y_continuous(expand = c(0,1), limits = c(0, RoundUpToNextZeros(max.frequency))) +
    theme_bw() +
    theme(plot.title = element_text(size = PLOT.TITLE.SIZE),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.title = element_text(size = AXIS.TITLE.SIZE),
          axis.text = element_text(size = AXIS.TEXT.SIZE))
  print(plot.quality.score)
}

dev.off()

