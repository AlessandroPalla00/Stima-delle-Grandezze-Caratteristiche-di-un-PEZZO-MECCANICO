# Stima delle Misure Caratteristiche di un Pezzo Meccanico tramite Triangolazione Laser 2D


Il progetto mira a caratterizzare dimensionalmente una campana dentata meccanica, stimandone i diametri caratteristici e valutando rigorosamente l'incertezza di misura attraverso simulazioni Monte Carlo, in accordo con le direttive GUM (Guide to the expression of Uncertainty in Measurement).

## Obiettivi del Progetto

L'obiettivo principale è la misurazione di tre grandezze geometriche critiche della campana dentata:

1.  D_in: Diametro interno alla radice dei denti (valli).
2.  D_out: Diametro esterno alla cresta dei denti.
3.  D_0: Diametro della circonferenza esterna, tangente a sfere ideali di raggio 4.25 mm incastrate nei fianchi dei denti.


## Architettura del Sistema di Misura (Setup Fisico)

L'acquisizione dei dati è stata effettuata tramite un banco di misura opto-meccanico composto da:

*   **Sensore Laser Triangolare 2D (Wenglor MLWL132):** Posizionato ortogonalmente all'asse di rotazione, acquisisce porzioni del profilo del pezzo.
*   **Rotatore Motorizzato (Physik Instrumente M-110.C8631):** Controllato via PC, permette la rotazione del pezzo a step angolari noti (30 gradi)
*   **Alimentatore Lineare a basso rumore (GW Instek GPS-3030D).**
*   **Supporto di centraggio:** Cilindro stampato in 3D per l'accoppiamento tra il rotatore e il pezzo in esame.

La strategia adottata è l'acquisizione a Campioni Finiti seguita da un Post-Processing in MATLAB. Le 12 scansioni parziali del profilo vengono allineate e fuse in un unico sistema di riferimento globale mediante trasformazioni roto-traslatorie, che tengono conto della discrepanza tra l'origine del laser e l'asse di rotazione meccanico.

## Analisi Metrologica e Modello delle Incertezze

Per validare il sistema di misura, è stato sviluppato un modello delle incertezze, classificando le fonti in due categorie:

1.  **Incertezza di Tipo A (Meccanica del Banco):**
    Deriva dai difetti cinematici del rotatore, in particolare l'eccentricità e lo sfarfallamento dell'asse durante la rotazione. È stata valutata acquisendo 13 pose di un cilindro di calibrazione. I centri geometrici di queste pose sono stati calcolati mediante il fit ai minimi quadrati di Kåsa. La deviazione standard della popolazione di questi centri (u_xc e u_yc) costituisce il parametro di dispersione per la distribuzione Gaussiana usata nella simulazione.

2.  **Incertezza di Tipo B (Limiti del Sensore):**
    Rappresenta i limiti intrinseci del sensore laser, estratti dal datasheet del costruttore:
    *   Linearità: 0.0325 mm
    *   Risoluzione: 0.014 mm

L'inclinazione residua del pezzo (tilt cosine error) dovuta al posizionamento manuale e alla messa a bolla rappresenta la principale incertezza di modello (Model Uncertainty), in quanto l'algoritmo sviluppato opera su base strettamente bidimensionale (2D).

## Metodologia di Calcolo e Motore Monte Carlo

Data la forte non-linearità degli algoritmi geometrici utilizzati, la propagazione analitica dell'incertezza (GUM standard) risulta inapplicabile. Si è pertanto adottato il Metodo Monte Carlo con 30.000 iterazioni eseguite in parallelo (`parfor`).

### Flusso dell'Algoritmo

1.  **Perturbazione del Setup:** Ad ogni iterazione, le coordinate del centro di rotazione nominale vengono perturbate utilizzando generatori casuali Gaussiani (Tipo A).
2.  **Roto-traslazione e Fusione:** Le 12 scansioni angolari (ruotate in fase di pre-caricamento per ottimizzare i tempi) vengono allineate e fuse nel nuovo sistema di riferimento perturbato.
3.  **Filtraggio e Fit di Kåsa (D_in e D_out):** La nuvola di punti completa viene suddivisa in 72 bin angolari. Attraverso l'ordinamento vettoriale nativo di MATLAB, vengono estratti i quantili per isolare le valli (root points) e le creste (tip points). Su questi set di punti viene applicato il Fit di Kåsa per determinare i diametri.
4.  **Estrazione Rette e Calcolo Tangenze (D_0):** I punti corrispondenti ai fianchi dei denti vengono isolati escludendo le zone di raccordo. Le rette dei fianchi vengono calcolate mediante l'analisi agli autovalori (PCA - eig). Successivamente, viene simulato l'incastro geometrico di una sfera da 8.5 mm tra le rette consecutive. I centri di queste sfere vengono interpolati con Kåsa per determinare D_0. L'operazione trigonometrica inversa per trovare l'incastro della sfera genera una forte asimmetria nella propagazione dell'errore.
5.  **Iniezione Errore Sensore:** Sui diametri finali calcolati viene iniettato l'errore Rettangolare del sensore laser (Tipo B).
6.  **Taglio Percentili Empirici:** Per evitare assunzioni sulla forma della distribuzione finale (spesso non Normale a causa delle non-linearità geometriche), l'Incertezza Estesa (U) al 95% viene estratta tagliando i percentili empirici al 2.5% (LCL) e 97.5% (UCL) sulla popolazione dei 30.000 risultati ordinati.

## Struttura del Codice

Il codice principale, `main_finale.m`, è strutturato nelle seguenti sezioni:

1.  **Parametri di Calibrazione:** Definizione del centro nominale e delle incertezze (Tipo A) per il setup corrente.
2.  **Setup Monte Carlo:** Generazione dei vettori perturbati (M=30.000).
3.  **Pre-caricamento Dati:** Importazione dei `.txt` parziali e pre-calcolo delle rotazioni.
4.  **Motore Monte Carlo Parallelo (`parfor`):** Esecuzione delle simulazioni con vettorializzazione avanzata (ordinamento `sort` al posto di funzioni pesanti come `quantile`) e calcolo rette tramite matrice di covarianza e `eig` in sostituzione della funzione built-in `pca`.
5.  **Aggiunta Incertezza Sensore:** Iniezione del Tipo B (Distribuzione Rettangolare).
6.  **Statistica e Taglio:** Estrazione media e limiti LCL/UCL empirici.
7.  **Plotting:** Stampa dei risultati e generazione dei grafici visivi (Profilo ricostruito e Istogrammi di probabilità).