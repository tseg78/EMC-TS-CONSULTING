EMC TS CONSULTING — Déclinaison BCI EMCTS V1 — application installable (PWA)

CONTENU
- index.html            : interface de l'application
- engine.js             : moteur de calcul (ligne multiconductrice, résolution nodale, auto-tests)
- manifest.webmanifest  : configuration d'installation
- service-worker.js     : mise en cache hors ligne
- icons/                : icônes 192 et 512 px (standard et « maskable »)

PUBLICATION SUR GITHUB PAGES
1. Placer tous les fichiers à la racine du dépôt (ou d'un dossier /docs), en conservant le dossier icons/.
2. Dans le dépôt : Settings > Pages > Source « Deploy from a branch », branche main, dossier / (root) ou /docs.
3. Ouvrir https://<utilisateur>.github.io/<dépôt>/ dans Chrome ou Edge, puis « Installer l'application »
   (bouton dans la page ou icône d'installation dans la barre d'adresse).
4. Sur iPhone : Safari > Partager > Sur l'écran d'accueil.
5. Après la première ouverture complète, l'application fonctionne hors ligne. Une nouvelle version
   publiée est prise en compte à l'ouverture suivante, dès que le poste est en ligne.

ESSAI EN LOCAL
Dans le dossier : python -m http.server 8000, puis http://localhost:8000/index.html
Un double-clic sur index.html permet aussi de tester l'interface, sans installation ni mode hors ligne.

CONFIGURATIONS
La configuration courante est conservée dans le navigateur. Les boutons « Exporter » et « Importer »
échangent des fichiers JSON pour archiver ou rejouer un cas d'essai.

STATUT DU MODÈLE
Modèle EMC TS CONSULTING en cours de validation. Le bouton « Lancer les vérifications du moteur »
(fenêtre Résultat) exécute onze contrôles analytiques. Les résultats doivent être confrontés à des mesures
avant toute utilisation comme exigence. Hypothèses détaillées dans la fenêtre Résultat.

DEUX MODES DE CALCUL
- Approché (par défaut) : le courant BCI total se répartit entre les fils selon les charges amont et aval
  de chacun. Harnais (5 cm du plan d'essai) et câblage interne (hauteur conducteur–châssis saisie) sont des
  lignes avec leur capacité parasite : les résonances apparaissent. Couplages entre fils négligés (majorant).
  Résistance linéique des conducteurs (cuivre, effet de peau) prise en compte par défaut dans les deux modes ;
  elle peut être négligée dans la fenêtre Harnais.
- Détaillé : lignes multiconductrices (géométrie du toron, position de la pince, pince calibrée), pour
  examiner les résonances du faisceau.

UTILISATION
À chaque ouverture, les blocs du synoptique sont blancs : les valeurs saisies sont conservées, mais chaque bloc
doit être relu puis appliqué. Un bloc passe au vert après « Appliquer » si ses paramètres sont complets.
Nombre de conducteurs : 60 au plus en mode approché, 40 en mode détaillé. Pour saisir vite de nombreux fils
identiques : les fils ajoutés reprennent la terminaison du dernier fil signal ou alimentation, et le bouton
« Copier vers tous les fils … » recopie une terminaison, un filtre ou une charge sur les fils de même fonction.
Le câblage interne est optionnel : absent par défaut (connecteur relié directement à la carte), il n'entre dans
le calcul que s'il est déclaré présent. Le résultat est disponible quand tous les blocs requis sont verts.

FAISCEAU BLINDÉ (modèle forfaitaire)
Blindage global « Blindé : atténuation forfaitaire » : le calcul est fait comme sans blindage, puis les courants et
tensions des fils intérieurs sont réduits de l'atténuation saisie (dB), identique à toutes les fréquences.
Approximation pour un blindage dont on ne connaît rien ; la f.é.m. requise n'est alors pas calculée.

CHARGES RÉSISTIVES RÉELLES
Pour une terminaison ou une charge carte de type R, deux champs facultatifs décrivent la résistance réelle :
inductance parasite série (nH) et capacité parasite parallèle (pF). Z = (R + jωLs) // Cp. Case vide = non modélisée.

EXPORT DES COURBES (fenêtre Résultat)
« Exporter les courbes (CSV) » : fréquence, courant total, f.é.m. requise, puis pour chaque fil le courant et la
tension sur la charge carte (V et dBµV). Format tableur français : séparateur « ; », virgule décimale, UTF-8.

MENU FICHIER (en haut à droite)
- Exporter : enregistre toute la configuration (niveaux, modèle, fils, filtres, charges) dans un fichier JSON,
  pour archiver un cas, le transmettre ou le rejouer plus tard.
- Importer : recharge une configuration exportée ; les blocs sont alors à réappliquer.
- Revenir au cas de référence : recharge l'exemple intégré à 3 fils (Signal 100 Ω, GND, Alimentation 10 Ω).

CHANGEMENTS PRINCIPAUX PAR RAPPORT AU PROTOTYPE V8.11
- Faisceau et câblage interne modélisés en lignes multiconductrices (L′ et C′ géométriques, mutuelles incluses),
  valables de 10 kHz à 400 MHz ; C′ cohérente avec la géométrie (plus de 470 pF/m imposés).
- Pince modélisée en source série au point d'injection saisi ; deux modes : courant imposé ou pince calibrée.
- Tous les champs de saisie pilotent le calcul (courant, filtres connecteur et carte, liaison 0 V carte,
  parasites HF des composants). Champs sans effet retirés.
- Faisceau blindé pris en compte par une atténuation forfaitaire.
- Bandes DO-160 Section 20 (10 kHz–400 MHz) et MIL-STD-461 CS114 (10 kHz–200 MHz).
- Export / import JSON, impression des résultats, icônes fournies.
