#import "@preview/amsterdammetje-article:0.1.1": abstract, article, heading-author

#set text(lang: "en")

#show: article(
  title: "Master thesis proposal",
  authors: "Wessel Beumer",
  ids: "12640662",
  link-outline: true,
  date: datetime(year: 2026, month: 04, day: 14),
)

= Base Data
+ Title: Expansion of MCM finding algorithms
+ Examiner: Alberto Pérez de Alba Ortíz
+ Assesor: Later
+ Daily supervisor: Clelia de Mulatier
+ Internal
+ Institute: ITFA

= Project description
Dr. Mulatier's group researches how spin models, such as the Ising model, can be used to model arbitrary data. This makes it possible to model and reveal high-order interactions in data in a more interpertable way than neural networks. They do this by splitting large models in collections of smaller models, called minimally complex models. The problem is that the programs they have written to find these models are limited to datasets with 128 variables. I am going to develop a new program which does not suffer from this constraint to enable the analysis of larger systems.

The extension to larger datasets has two main challenges. (1) The efficiency of the code relies on encoding data and models over 128-bit integers to perform bit operations (hence the limitations to 128 variables). To extend to larger systems, this data structure must change; we aim to find a solution that would still retain the efficiency of using bit operations (or almost). (2) The space of possible models increases superexponentially with the number of variables, making the optimization problem more challenging. This means that we need to find more efficient algorithms if we want to still find near-optimal models. We will explore different algorithms (such as greedy- or simulated annealing-based, and parallel tempering). We may also look at implementing a Potts model version of this technique, depending on the time available.

= Timetable
+ Start date: 2026-03-30
+ Planned end date: 2026-11-30
+ Interuptions: Im on vacation 9-14 of august
