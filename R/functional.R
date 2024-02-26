# =====functional gene analysis====
{
    COG_desc <- data.frame(
        Category = c(
            "INFORMATION STORAGE AND PROCESSING", "INFORMATION STORAGE AND PROCESSING", "INFORMATION STORAGE AND PROCESSING", "INFORMATION STORAGE AND PROCESSING", "INFORMATION STORAGE AND PROCESSING",
            "CELLULAR PROCESSES AND SIGNALING", "CELLULAR PROCESSES AND SIGNALING", "CELLULAR PROCESSES AND SIGNALING", "CELLULAR PROCESSES AND SIGNALING", "CELLULAR PROCESSES AND SIGNALING",
            "CELLULAR PROCESSES AND SIGNALING", "CELLULAR PROCESSES AND SIGNALING", "CELLULAR PROCESSES AND SIGNALING", "CELLULAR PROCESSES AND SIGNALING", "CELLULAR PROCESSES AND SIGNALING",
            "CELLULAR PROCESSES AND SIGNALING", "METABOLISM", "METABOLISM", "METABOLISM", "METABOLISM", "METABOLISM", "METABOLISM", "METABOLISM", "METABOLISM",
            "POORLY CHARACTERIZED", "POORLY CHARACTERIZED"
        ),
        Code = c("J", "A", "K", "L", "B", "D", "Y", "V", "T", "M", "N", "Z", "W", "U", "O", "X", "C", "G", "E", "F", "H", "I", "P", "Q", "R", "S"),
        Description = c(
            "Translation, ribosomal structure and biogenesis", "RNA processing and modification", "Transcription", "Replication, recombination and repair", "Chromatin structure and dynamics",
            "Cell cycle control, cell division, chromosome partitioning", "Nuclear structure", "Defense mechanisms", "Signal transduction mechanisms", "Cell wall/membrane/envelope biogenesis",
            "Cell motility", "Cytoskeleton", "Extracellular structures", "Intracellular trafficking, secretion, and vesicular transport", "Posttranslational modification, protein turnover, chaperones",
            "Mobilome: prophages, transposons", "Energy production and conversion", "Carbohydrate transport and metabolism", "Amino acid transport and metabolism", "Nucleotide transport and metabolism",
            "Coenzyme transport and metabolism", "Lipid transport and metabolism", "Inorganic ion transport and metabolism", "Secondary metabolites biosynthesis, transport and catabolism",
            "General function prediction only", "Function unknown"
        ),
        stringsAsFactors = FALSE
    )
}
