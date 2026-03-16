#NLES5 Functions 
# L function NLES5 -----
nles5 <- function(Y, #Temporal tendency
                  Ntheta, # Nitrogen
                  C, # Crop
                  P, # Perc
                  # Psas, # Perc sas code
                  S, # Soil
                  # EEA=0, Fdato=1, EMA=0.05, ETS=0.05, EPJ=1/11, # NUAR coef
                  tao = -0.1108,# Parameter for time effect
                  mu = 23.51,# Parameter for base level se =4.3418
                  k = 1.5,# Parameter for scaling C and N
                  rho = 1.085, # Scaling factor accounting for bias transformation this is a problem
                  varres= 3.6274910633 # from sas25/08
) {
  # Calculate the nitrogen leaching using the N function with provided parameters
  #𝐿 = 𝜏(𝑌 – 1991) + {(𝜇 + 𝜃𝑖𝑁 + 𝐶)^𝜅}(𝑃 𝑆)p
  
  L <-  tao * (Y - 1991) + ((mu + Ntheta + C)^k)*((P * S) * rho)
  
  # Calculate L using the percolation function as found in SAS code
  # L_percSAS <-  tao * (Y - 1991) + ((mu + Ntheta + C)^k)*(Psas * S) * rho
  # Calculate L by rho affecting entire expression as in SAS code
  Lwr <- (tao * (Y - 1991) + sqrt(((mu + Ntheta + C)^k)*(P * S)))^2+varres #(tao * (Y - 1991) + ((mu + Ntheta + C)^k)*(P * S)) * rho
  
  #𝐿NUAR = (𝜏(𝑌 − 1991) + {(𝜇 + 𝜃𝑖𝑁 + 𝐶)𝜅}(𝑃 𝑆)p)(1 - EEA Fdato - EMA - ETS)(1- EPJ)
  # L_nuar <- L * (1 - EEA*Fdato - EMA - ETS) * (1 - EPJ)
  
  return(list(L=L,
              Lwr=Lwr,
              #L_nuar=L_nuar,
              Ntheta=Ntheta,
              C=C,
              P=P,
              #Psas=Psas,
              S=S#,
              
              #L_percSAS=L_percSAS
  ))
  
}

# N component ----
# 1.1 Define the Nitrogen Model Function
# The model calculates N based on various nitrogen parameters and input variables.
# N = βt NT + βCS MNCS + βCA MNCA + βudb MNudb + βm1 (M1+M2)/2 + βf0 F0 +
#     βf1 (F1+F2)/2 + βg0 G0 + βm1 (G1+G2)/2

N_func <- function(
    #Inputs
  NT, # Total N in topsoil (0-25 cm) (ton N)
  MNCS, # Mineral N spring in harvest year
  MNCA, # Mineral N autumn in harvest year
  MNudb, # Mineral N from grazing in harvest year
  M1, M2, # N mineral prev years
  F0, F1, F2, #N from N fixation
  G0, G1, G2, #N from Organic fertilization
  WC,  # autum vegetation
  
  #define param
  beta_t = 0.456793,    # Total N in topsoil (0-25 cm) (ton N)
  beta_CS = 0.049570,    # Mineral N spring in harvest year
  beta_CA = 0.157044,   # Mineral N autumn in harvest year
  beta_udb = 0.038245,  # Mineral N from grazing in harvest year
  beta_m1_M = 0.026499, # Mineral N added to previous/pre-previous crop (M1+M2)
  beta_f0 = 0.016314,   # Nitrogen fixation in harvest year
  beta_f1 = 0.026499,    # Nitrogen fixation in previous/pre-previous crop (F1+F2)
  beta_g0 = 0.014099,   # Organic N spring in harvest year
  beta_m1_G = 0.026499,   # Organic N added to previous/pre-previous crop (G1+G2)
  theta_2 = 1.205144 # Factor for N in winter crops (WC=2) or not (WC=1), default is 1.205144
  
) {
  
  
  # Calculate the N value based on the provided parameters and input variables
  N <- (
    beta_t * NT +
      beta_CS * MNCS +
      beta_CA * MNCA +
      beta_udb * MNudb +
      beta_m1_M * ((M1 + M2) / 2) +
      beta_f0 * F0 +
      beta_f1 * ((F1 + F2) / 2) +
      beta_g0 * G0 +
      beta_m1_G * ((G1 + G2) / 2)
  )
  Ntheta <- ifelse(WC==1,N, N * theta_2)
  return(Ntheta)
}

# C component ----


C_func <- function(M, W, MP, WP) {
  # Parameter tables
  M_p <- c(
    M1 = 0,              # Vintersæd (Winter cereals)
    M2 = -6.744,         # Vårsæd (Spring cereals)
    M3 = -7.279,         # Bælgsæd-korn blanding (Legume-cereal mix)
    M4 = -13.493,        # Græs og kløvergræs (Grass and clover grass)
    M5 = -17.478,        # Frøgræs (Seed grass)
    M6 = -11.192,        # Brak (Fallow)
    M8 = -0.640,         # Sukkerroer og foderroer (Sugar beets and fodder beets)
    M9 = 3.534,          # Majshelsæd og kartofler (Maize silage and potatoes)
    M10 = -7.319,         # Vinterraps (Winter rapeseed)
    M11 = -1.248,        # Vintersæd efter græs, kløvergræs, frøgræs og brak (Winter cereals after grass, clover, seed grass, or fallow)
    M12 = 19.524,        # Majshelsæd efter græs, kløvergræs, frøgræs og brak (Maize silage after grass, clover, seed grass, or fallow)
    M13 = -6.229,        # Vårsæd efter græs, kløvergræs, frøgræs og brak (Spring cereals after grass, clover, seed grass, or fallow)
    M14 = -2.866         # Bælgsæd og vårraps (Legumes and spring rapeseed)
  )
  
  W_p <- c(
    W1 = 0,              # Vintersæd (Winter cereals)
    W2 = -2.055,         # Bar jord (Bare soil)
    W3 = -0.456,         # Bar jord efter majshelsæd og kartofler (Bare soil after maize silage and potatoes)
    W4 = -15.959,        # Efterafgrøder, undersået græs og brak (Catch crops, undersown grass, and fallow)
    W5 = -3.792,         # Ukrudt og spildkornsplanter (Weeds and volunteer cereals)
    W6 = -14.596,        # Græs, kløvergræs, vinterraps, roer, frøgræs og udlægsafgrøder (Grass, clover, winter rapeseed, beets, seed grass, and undersown crops)
    W8 = -21.060,         # Vintersæd efter græs og frøgræs (Winter cereals after grass and seed grass)
    W9 = -1.049         # Græs og kløvergræs pløjet sent efterår eller vinter (Grass and clover plowed late autumn or winter)
  )
  
  MP_p <- c(
    MP1 = 0,             # Vintersæd (Winter cereals)
    MP2 = 2.847,         # Andre afgrøder end vintersæd, græs, kløvergræs, frøgræs og brak (Other crops)
    MP3 = 0.664,         # Græs, kløvergræs, frøgræs og brak (Grass, clover, seed grass, and fallow)
    MP4 = 1.160,          # Vår- og vinterafgrøder efter græs, kløvergræs, frøgræs og brak (Spring/winter crops after grass, clover, seed grass, or fallow)
    MP12 = 2.847
  )
  
  WP_p <- c(
    WP1 = 0,             # Vintersæd (Winter cereals)
    WP2 = 9.704,         # Bar jord og spildkorn (Bare soil and volunteer cereals)
    WP3 = 10.601,        # Græs og kløvergræs (Grass and clover)
    WP4 = 9.354,         # Efterafgrøder (Catch crops)
    WP5 = 13.241,        # Frøgræs og brak (Seed grass and fallow)
    WP6 = 5.483,         # Sukkerroer, foderroer og hamp (Sugar beets, fodder beets, and hemp)
    WP7 = -1.572,        # Bar jord efter majshelsæd og kartofler (Bare soil after maize silage and potatoes)
    WP8 = 7.413,         # Vinterraps (Winter rapeseed)
    WP9 = 7.396,         # Bar jord eller vintersæd efter græs, kløvergræs, frøgræs eller brak ompløjet *forår* (Bare soil or winter cereals after spring-plowed grass, clover, seed grass, or fallow)
    WP10 = 10.975        # Bar jord eller vintersæd efter græs, kløvergræs, frøgræs eller brak ompløjet *efterår* (Bare soil or winter cereals after autumn-plowed grass, clover, seed grass, or fallow)
  )
  
  #Look up values
  C <- M_p[[paste("M", M, sep = "")]] + #paste("M",M[[1]], sep = "")
    W_p[[paste("W", W, sep = "")]] +
    
    MP_p[[paste("MP", MP, sep = "")]] +
    WP_p[[paste("WP", WP, sep = "")]]
  
  return(C)
}

# S component ----
S_func <-  function(P_ler=0.001849, CU){
  
  S <- exp(P_ler*CU)
  
  return(S)
}

# P component documentation ----
P_func <- function(jbnr,
                   AAa, #d1 April–August current year
                   AAb, #d2 Sept–March current year
                   APb, #d3   Sept–March preceding year
                   
                   # Parameters for sandy soil (JB1 + JB3)
                   delta1s = 0.001194,     # April–August percolation on sandy soil
                   delta2s = 0.001107,     # Sept–March percolation on sandy soil
                   noo2s = 0.000856,  # Sept–March percolation in preceding year on sandy soil
                   
                   # Parameters for loamy soil (all other JB types)
                   delta1c = 0.000798,     # April–August percolation on loamy soil
                   delta2c = 0.000745,     # Sept–March percolation on loamy soil
                   noo2c = 0.000638  # Sept–March percolation in preceding year on loamy soil
) {
  
  # Compute P based on soil type
  if (jbnr <= 3) {
    # Sandy soil formula
    P <- (1 - exp(-delta1s * AAa - delta2s * AAb)) * exp(-noo2s * APb)
  } else {
    # Loamy soil formula
    P <- (1 - exp(-delta1c * AAa - delta2c * AAb)) * exp(-noo2c * APb)
  }
  
  return(P)
}

# P component SAS script ----
Psas_func <- function(jbnr,
                      AAa, #d1 April–August current year
                      AAb, #d2 Sept–March current year
                      APb, #d3   Sept–March preceding year
                      p2,
                      p3,
                      
                      # Parameters for sandy soil (JB1 + JB3)
                      delta1s = 0.001194,     # April–August percolation on sandy soil
                      delta2s = 0.001107,     # Sept–March percolation on sandy soil
                      noo2s = 0.000856,  # Sept–March percolation in preceding year on sandy soil
                      
                      # Parameters for loamy soil (all other JB types)
                      delta1c = 0.000798,     # April–August percolation on loamy soil
                      delta2c = 0.000745,     # Sept–March percolation on loamy soil
                      noo2c = 0.000638  # Sept–March percolation in preceding year on loamy soil
) {
  
  if (jbnr <= 3) {
    # Sandy soil formula for Psas
    Psas <- (1 - exp(-delta1s * AAa - delta2s * AAb - delta2s * APb)) * exp(-noo2s * p2 - noo2s * p3)
  } else {
    # Loamy soil formula for Psas
    Psas <-  (1 - exp(-delta1c * AAa - delta2c * AAb - delta2c * APb)) * exp(-noo2c * p2 - noo2c * p3)
  }
  
  return(Psas)
}

# Percolation exercise ==========================================
library(tidyverse)
library(readxl)
library(ggrepel)
 
# Marginal effects NLES5: -----
# We vary each input variable across its observed range while keeping all other variables fixed at their mean (for continuous) or mode (for categorical). This allows us to isolate the effect of each variable on the predicted leaching (L) and visualize it in a faceted plot.

## --- 1. Load and Prepare Clean Base Data ----
# We map the raw dataset directly to the variable names expected by the functions 
# so we don't need complex translation vectors later.
sasoutput <- read_excel("Scenarier20190909B4_found0325.xls") |>
  rename("na"="NA")

base_data <- sasoutput |>
  transmute(
    Y = Indb_aar,
    
    # N components
    NT = TN, MNCS = NS, MNCA = na, MNudb = Nudb,
    M1 = nlevelMin1, M2 = nlevelMin2, F0 = NlevelFix0, F1 = NlevelFix1, F2 = NlevelFix2,
    G0 = NlevelGod0, G1 = NlevelGod1, G2 = NlevelGod2, WC = Vafgr_Kappa,
    
    # C components
    M = Mau, MP = Mfu, W = Vau, WP = Vfu,
    
    # P components
    AAa = d1, 
    AAb = d2 + d3, 
    APb = p2 + p3, 
    jbnr = jbnr,
    
    # S components
    CU = CU
  ) |> drop_na()

## --- 2. Helper functions & Fixed Values ----
get_mode <- function(v) {
  uniqv <- unique(na.omit(v))
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# Define which variables are categorical
categorical_vars <- c("M", "W", "MP", "WP", "WC", "jbnr", "SoilG")

# Calculate mean (for continuous) or mode (for categorical) to act as our "fixed" baseline
fixed_values <- list()
for (col_name in names(base_data)) {
  if (col_name %in% categorical_vars) {
    fixed_values[[col_name]] <- get_mode(base_data[[col_name]])
  } else {
    fixed_values[[col_name]] <- mean(base_data[[col_name]], na.rm = TRUE)
  }
}

## --- 3. Component and Color Definitions ----
component_mapping <- c(
  "Y"="Year",
  "MNCS" = "Nitrogen", "MNCA" = "Nitrogen", "M1" = "Nitrogen", "M2" = "Nitrogen",
  "F0" = "Nitrogen", "F1" = "Nitrogen", "F2" = "Nitrogen", "G0" = "Nitrogen",
  "G1" = "Nitrogen", "G2" = "Nitrogen", "MNudb" = "Nitrogen", "NT" = "Nitrogen", "WC" = "Nitrogen",
  "M" = "Crop", "MP" = "Crop", "W" = "Crop", "WP" = "Crop",
  "AAa" = "Percolation", "AAb" = "Percolation", "APb" = "Percolation", "jbnr" = "Percolation",
  "CU" = "Soil"
)

component_colors <- c(
  "Year" = "#FFAE00",
  "Soil" = "#e45a3d",
  "Nitrogen" = "#377eb8",
  "Crop" = "#3AA600",
  "Percolation" = "#985aa1"
)

## --- 4. Perform Marginal Effect Analysis ----
# run functions in 1_fgk_NLES5.R to environment first so we can call them in the loop without overhead
marginal_effects_results <- list()
input_variables_to_test <- names(base_data)

for (var_to_vary in input_variables_to_test) {
  
  # Create a sequence of values to test for the current variable
  if (var_to_vary %in% categorical_vars) {
    # Test all unique classes for categorical variables
    test_vals <- sort(unique(base_data[[var_to_vary]]))
  } else {
    # Test 100 evenly spaced points from min to max for continuous variables
    test_vals <- seq(min(base_data[[var_to_vary]]), max(base_data[[var_to_vary]]), length.out = 100)
  }
  
  # Build a dataframe where everything is fixed EXCEPT the variable we are testing
  temp_df <- as.data.frame(lapply(fixed_values, function(x) rep(x, length(test_vals))))
  temp_df[[var_to_vary]] <- test_vals
  
  # Run the model using dplyr::rowwise() (Much cleaner than a nested for loop)
  results_df <- temp_df |>
    rowwise() |>
    mutate(
      N_est = N_func(NT=NT, MNCS=MNCS, MNCA=MNCA, MNudb=MNudb, M1=M1, M2=M2, F0=F0, F1=F1, F2=F2, G0=G0, G1=G1, G2=G2, WC=WC),
      C_est = C_func(M=M, W=W, MP=MP, WP=WP),
      P_est = P_func(jbnr=jbnr, AAa=AAa, AAb=AAb, APb=APb),
      S_est = S_func(CU=CU),
      L = nles5(Y=Y, Ntheta=N_est, C=C_est, P=P_est, S=S_est)$L,
      
      # Tagging metadata for plotting
      varied_variable = var_to_vary,
      original_value = .data[[var_to_vary]],
      Component = component_mapping[var_to_vary],
      Is_Categorical = var_to_vary %in% categorical_vars
    ) |>
    ungroup() |>
    select(varied_variable, original_value, L, Component, Is_Categorical)
  
  marginal_effects_results[[var_to_vary]] <- results_df
}

# Combine all results
final_results_df <- bind_rows(marginal_effects_results) |>
  mutate(
    varied_variable = factor(varied_variable, levels = names(component_mapping)),
    Component = factor(Component, levels = names(component_colors))
  )

## --- 5. Plotting: Faceted Marginal Effects ----

# We split the plot logic slightly: Lines for continuous, Points for categorical
p_marginal <- ggplot(final_results_df, 
                     aes(x = original_value, y = L, color = Component)) +
  # Add lines for continuous variables
  geom_line(data = filter(final_results_df, !Is_Categorical), linewidth = 1) +
  # Add points for categorical variables
  geom_point(data = filter(final_results_df, Is_Categorical), size = 2.5) +
  labs(
    title = "Marginal Effects of NLES5 Parameters",
    subtitle = "Showing the effect on Annual Leaching (L) while all other inputs are fixed at their mean/mode",
    x = "Input Variable Value",
    y = "Predicted Annual Leaching (L)"
  ) +
  facet_wrap(~ varied_variable, scales = "free_x", ncol = 5) +
  scale_color_manual(values = component_colors) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 9, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey40"),
    legend.position = "bottom"
  )

print(p_marginal)
# 
# ggsave("marginal_effects_nles5.png", plot = p_marginal, 
#        units = "mm", width = 250, height = 200, dpi = 300)


## --- 5. Plotting: Separated Percolation Marginal Effects ----

### 5.1 Continuous Variables Plot ----
p_marginal_cont <- final_results_df |>
  filter(!Is_Categorical & Component=="Percolation") |> 
  ggplot(aes(x = original_value, y = L, color = Component)) +
  geom_line(linewidth = 1) +
  labs(
    title = "Marginal Effects: Continuous Parameters",
    subtitle = "Effect on Annual Leaching (L) while other inputs are fixed at mean/mode",
    x = "Input Variable Value",
    y = "Predicted Annual Leaching (L)"
  ) +
  facet_wrap(~ varied_variable, scales = "free_x", ncol = 4) +
  scale_color_manual(values = component_colors) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey40"),
    legend.position = "bottom"
  )

print(p_marginal_cont)
ggsave("marginal_effects_continuous.png", plot = p_marginal_cont, units = "mm", width = 250, height = 200, dpi = 300)


### 5.2 Categorical Variables Plot ----

# Map the numeric codes to crop names
recode_map <- c(
  "M-1" = "Winter cereal", "M-2" = "Spring cereal", "M-3" = "Grain-legume mix", 
  "M-4" = "Grass and grass-clover", "M-5" = "Grass for seed", "M-6" = "Set-aside", 
  "M-8" = "Beet and hemp", "M-9" = "Maize and potato", "M-10" = "Winter Oilseed rape", 
  "M-11" = "Winter cereal after grass", "M-12" = "Maize after grass", 
  "M-13" = "Spring cereal after grass", "M-14" = "Grain legume & Spring oilseed rape",
  
  "MP-1" = "Winter cereal", "MP-2" = "Other crops", "MP-3" = "Grass or grass-clover", 
  "MP-12" = "Spring or Winter crops after grass",
  
  "W-1" = "Winter cereal", "W-2" = "Bare Soil", "W-3" = "Autumn cultivation", 
  "W-4" = "Cover crop", "W-5" = "Weeds, volunteers", "W-6" = "Winter oilseed rape, grass-clover", 
  "W-8" = "Winter cereal after grass", "W-9" = "Grass ploughed late autumn/winter",
  
  "WP-1" = "Winter cereal", "WP-2" = "Bare soil", 
  "WP-3" = "Grass-clover and Winter oilseed rape ", "WP-4" = "Cover grass", 
  "WP-5" = "Grass for seeds and set aside", "WP-6" = "Beets and hemp", 
  "WP-7" = "Bare soil after maize and potato", 
  "WP-11" = "Grass followed by spring cereal, maize or beets",
  
  "WC-1" = "Crops w/ large N uptake", "WC-2" = "Crops w/ low/moderate N uptake",
  "jbnr-3" = "Sandy soil (JB 1-3)", "jbnr-4" = "Loamy soil (JB >3)"
)

# Prepare data and apply names
cat_plot_data <- final_results_df |>
  filter(Is_Categorical) |>
  mutate(
    # Create a unique key like "M-1" or "WP-11" to match the dictionary
    temp_key = paste0(varied_variable, "-", as.character(original_value)),
    Cat_Name = recode(temp_key, !!!recode_map, .default = as.character(original_value))
  )

p_marginal_cat <- cat_plot_data |>
  filter(Component=="Percolation") |> 
  #filter(Component=="Crop") |> 
  # reorder ensures the categories appear in their original numeric order rather than alphabetically
  ggplot(aes(x = reorder(Cat_Name, original_value), y = L, color = Component)) +
  geom_point(size = 3.5) +
  geom_segment(aes(x = reorder(Cat_Name, original_value), xend = reorder(Cat_Name, original_value), y = 0, yend = L), linetype = "dotted") +
  labs(
    title = "Marginal Effects: Categorical Parameters",
    subtitle = "Effect on Annual Leaching (L) while other inputs are fixed at mean/mode",
    x = "",
    y = "Predicted Annual Leaching (L)"
  ) +
  facet_wrap(~ varied_variable, scales = "free", ncol = 2) +
  scale_color_manual(values = component_colors) +
  coord_flip() + # Flips the plot so long crop names are perfectly readable on the Y axis
  theme_light() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 9),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey40"),
    legend.position = "bottom",
    panel.grid.major.y = element_blank() # Clean up horizontal lines since we use dotted segments
  )

print(p_marginal_cat)
# ggsave("marginal_effects_categorical.png", plot = p_marginal_cat, units = "mm", 
#        width = 250, height = 280, dpi = 300)

# Exploratory percolation distribution ----

# Define the requested Mustard (Sandy) and Reddish-Brown (Loamy) colors
soil_palette <- c(
  "Sandy (JB <= 3)" = "#E1AD01", 
  "Loamy (JB > 3)"  = "#8F4513"  
)

base_data |> 
  mutate(JB =ifelse(jbnr <= 3, "Sandy (JB <= 3)","Loamy (JB > 3)")) |> 
  pivot_longer(cols = c(AAa, AAb, APb), 
               names_to = "Percolation_Var", values_to = "Value") |> 
  ggplot(aes(x = Value, fill = JB)) +
  geom_density(position = "dodge", bins = 30, alpha = 0.7) +
  facet_wrap(~ Percolation_Var, scales = "free_x") +
  labs(
    title = "Distribution of Percolation Variables",
    x = "Percolation Value (mm)",
    y = "Frequency",
    fill = "Variable"
  ) +
  scale_fill_manual(values = soil_palette) +
  theme_minimal()

base_data |> 
  mutate(JB =ifelse(jbnr <= 3, "Sandy (JB <= 3)","Loamy (JB > 3)")) |> 
  pivot_longer(cols = c(AAa, AAb, APb), 
               names_to = "Percolation_Var", values_to = "Value") |> 
  ggplot(aes(y = Value, x=JB, fill = JB)) +
  geom_boxplot(position = "dodge", bins = 30, alpha = 0.7) +
  facet_wrap(~ Percolation_Var, scales = "free_x") +
  labs(
    title = "Distribution of Percolation Variables",
    y = "Percolation Value (mm)",
    x = "JB Category",
    fill = "JB Category"
  ) +
  scale_fill_manual(values = soil_palette) +
  theme_minimal()+
  theme(legend.position = "none")

# Marginal effects for percolation variables only, split by soil type (jbnr) ----
## --- 1. Filter for Percolation Continuous Variables ----
percolation_vars <- c("AAa", "AAb", "APb")
jbnr_levels <- c(3, 4) # 0 = Sandy (<=3), 1 = Loamy (>3)

perc_sensitivity_results <- list()

for (var_to_vary in percolation_vars) {
  
  # Create the test range for the percolation variable
  test_vals <- seq(min(base_data[[var_to_vary]]), 
                   max(base_data[[var_to_vary]]), 
                   length.out = 100)
  
  # Create a grid: every test value paired with BOTH jbnr levels
  # Everything else remains fixed at mean/mode
  temp_df <- as.data.frame(lapply(fixed_values, function(x) rep(x, length(test_vals) * 2)))
  
  # Overwrite the specific variables for this comparison
  temp_df[[var_to_vary]] <- rep(test_vals, each = 2)
  temp_df$jbnr <- rep(jbnr_levels, times = length(test_vals))
  
  # Run the model
  results_df <- temp_df |>
    rowwise() |>
    mutate(
      N_est = N_func(NT=NT, MNCS=MNCS, MNCA=MNCA, MNudb=MNudb, M1=M1, M2=M2, 
                     F0=F0, F1=F1, F2=F2, G0=G0, G1=G1, G2=G2, WC=WC),
      C_est = C_func(M=M, W=W, MP=MP, WP=WP),
      P_est = P_func(jbnr=jbnr, AAa=AAa, AAb=AAb, APb=APb),
      S_est = S_func(CU=CU),
      L = nles5(Y=Y, Ntheta=N_est, C=C_est, P=P_est, S=S_est)$L,
      
      varied_variable = var_to_vary,
      original_value = .data[[var_to_vary]],
      soil_type = ifelse(jbnr <= 3, "Sandy (JB <= 3)", "Loamy (JB > 3)")
    ) |>
    ungroup() |>
    select(varied_variable, original_value, L, P_est, soil_type)
  
  perc_sensitivity_results[[var_to_vary]] <- results_df
}

# Combine results
perc_comparison_df <- bind_rows(perc_sensitivity_results)

# --- 2. Plotting the Split Curves ---
p_perc_split <- perc_comparison_df |> 
  ggplot(aes(x = original_value, y = L, color = soil_type)) +
  geom_line(linewidth = 1.2) +
  labs(
    title = "Percolation Sensitivity split by Soil Type (jbnr)",
    subtitle = "Marginal effects on Leaching (L) for Sandy vs Loamy soils",
    x = "Percolation Value (mm)",
    y = "Predicted Annual Leaching (L)",
    color = "Soil Category"
  ) +
  facet_wrap(~ varied_variable, scales = "free_x") +
  scale_color_manual(values = c("Sandy (JB <= 3)" = "#D95F02", "Loamy (JB > 3)" = "#1B9E77")) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 11)
  )

print(p_perc_split)

## --- 2. Plotting the Split Curves ----
p_perc_split_P <- perc_comparison_df |> 
  ggplot(aes(x = original_value, y = P_est, color = soil_type)) +
  geom_line(linewidth = 1.2) +
  labs(
    title = "Percolation Sensitivity split by Soil Type (jbnr)",
    subtitle = "Marginal effects on Percolation term for Sandy vs Loamy soils",
    x = "Percolation Value (mm)",
    y = "P term Estimate",
    color = "Soil Category"
  ) +
  facet_wrap(~ varied_variable, scales = "free_x") +
  scale_color_manual(values = c("Sandy (JB <= 3)" = "#D95F02", "Loamy (JB > 3)" = "#1B9E77")) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 11)
  )

print(p_perc_split_P)

### Point wise check: Is the model actually sensitive to jbnr? ----
P_func(jbnr = 1, AAa = 400, AAb = 100, APb = 50)
P_func(jbnr = 4, AAa = 400, AAb = 100, APb = 50)


## 3. Gap Analysis ----
# Data processing for Gap and Slope
soil_gap <- perc_comparison_df |>
  select(varied_variable, original_value, L, soil_type) |>
  pivot_wider(names_from = soil_type, values_from = L) |>
  rename(L_sandy = `Sandy (JB <= 3)`, L_loamy = `Loamy (JB > 3)`) |>
  mutate(
    L_diff = L_sandy - L_loamy,
    slope_sandy = (L_sandy - lag(L_sandy)) / (original_value - lag(original_value)),
    slope_loamy = (L_loamy - lag(L_loamy)) / (original_value - lag(original_value))
  ) |>
  filter(!is.na(slope_sandy))

# Plot 2: The Gap (Penalty)
p_gap <- ggplot(soil_gap, aes(x = original_value)) +
  geom_ribbon(aes(ymin = 0, ymax = L_diff), fill = "#E1AD01", alpha = 0.3) +
  geom_line(aes(y = L_diff), color = "#8B4513", linewidth = 1) +
  facet_wrap(~ varied_variable, scales = "free_x") +
  labs(
    title = "The 'Sandy Soil Penalty'",
    subtitle = "Difference in Leaching (L_sandy - L_loamy) across percolation ranges",
    x = "Percolation Value (mm)",
    y = "Additional Leaching (kg N/ha)"
  ) +
  theme_minimal()

print(p_gap)

## 4.  Slopes Analysis  Compare the"jumps" (Sensitivities) ----
# This explains IF one soil type reacts more aggressively to rain than the other
slope_long <- soil_gap |>
  select(varied_variable, original_value, slope_sandy, slope_loamy) |>
  pivot_longer(cols = starts_with("slope"), names_to = "soil_type", values_to = "slope_val")

p_slope_comp <- ggplot(slope_long, aes(x = original_value, y = slope_val, color = soil_type)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ varied_variable, scales = "free") +
  scale_color_manual(
    values = c("slope_sandy" = "#E1AD01", "slope_loamy" = "#8B4513"),
    labels = c("Sandy Sensitivity", "Loamy Sensitivity")
  ) +
  labs(
    title = "Marginal Sensitivity Comparison",
    subtitle = "How many kg N/ha are gained/lost per 1mm of additional percolation?",
    x = "Percolation Value (mm)",
    y = "Slope (ΔL / ΔPerc)",
    color = "Legend"
  ) +
  theme_minimal() + theme(legend.position = "bottom")

print(p_slope_comp)

## 5 Maximum differences ----
# Final Summary Stats for Interpretation
summary_results <- soil_gap |>
  group_by(varied_variable) |>
  summarize(
    Avg_Sandy_Sens = mean(slope_sandy),
    Avg_Loamy_Sens = mean(slope_loamy),
    Max_Penalty = max(L_diff),
    Mean_Penalty = mean(L_diff)
  )

print(summary_results)
