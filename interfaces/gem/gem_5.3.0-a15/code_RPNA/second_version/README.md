# CodeApSwitch
Ici, on utilise le code original **gem_5.3.0-a15** ou **scm_2.3.0-a15** avec les modifications 
decrites dans le document:

[https://007gc-my.sharepoint.com/:w:/r/personal/frederick_chosson_ec_gc_ca/Documents/PC3_RPNA/De_CTL_a_Seamless.docx?d=w876db081dd974bd18e2f5bc70a7c541e&csf=1&web=1&e=Ugn1iH](De_CTL_a_Seamless)

Le but est de connaitre les impacts des 3 differentes parametrisations de autoconversion dans P3
La base commune est le "CodeA" sur la version a15 de GEM. Ce code NE PERMET PAS utilisation de SCPF.

Ici on utilise un Switch dans gem_settings.nml pour la parametrisation de autoconversion/accretion/self-collection : iparam=3 par defaut

On a rajouter l option supid pour le seuil de sursaturation p.r.a la glace avant d'appeler Cooper (approche de Paul Vaillancourt).

On trouve le code compile ici:
`~/home/gem/5.3.0-a15/CodeApSwitch`

Ce code est identique au code `~/LIENS/CodesP3_SCPF/original_123aFix` plus le Switch sur iparam

