Hallo Jörg und Tobias,

im Anhang schicke ich euch die neueste Version der RuNNer-LAMMPS-Implementierung, die wesentliche Geschwindigkeitsverbesserungen enthält. Installation und Verwendung sind unverändert (getestet mit LAMMPS 28Jun14), intern hat sich nur das Handling der Symmetriefunktion verändert. Diese sind jetzt in Gruppen zu je gleichen Elementkombinationen und Cutoff-Radien zusammengefasst, d.h. innerhalb jeder Gruppe sind die gleichen Nachbaratome relevant. Dadurch müssen viele Zwischenergebnisse (Abstände, Winkel,...) nicht für jede Symmetriefunktion erneut berechnet werden.

Bei den winkelabhängigen Symmetriefunktionen habe ich die kostspielige Powerfunktion pow() durch eine schnellere Integer-Version ersetzt, die automatisch verwendet wird, wenn der Exponent (zeta) eine Ganzzahl ist. Außerdem wird auch das Resultat der Exponentialfunktion exp(-eta*...) bei gleichem eta innerhalb jeder Gruppe wiederverwertet, d.h. zusätzliche Symmetriefunktionen mit gleichem eta sollten nicht mehr so viel Berechnungszeit kosten.

Ich habe die Beispielsysteme im "example"-Ordner und H2O ausprobiert und folgende Speedups erhalten:

H2O: 2.7x (single core) bis zu 3.8x (256 cores)
Cu: 9.5x (single core)
Cu2S: 5.3x (single core)
CuZnO: 7.0x (single core)

Die alte Variante der Symmetriefunktionsberechnung ist nach wie vor lauffähig enthalten. Man erhält sie, wenn man in Zeile 122 in "pair_runner.cpp"

calc_all_symfunc_group();

durch

calc_all_symfunc();

ersetzt. Bei meinen Tests lief alles problemlos, ich hoffe es sind keine Bugs mehr drinnen. Gebt mir bitte Bescheid falls ihr doch etwas findet.

Liebe Grüße,

Andi 
