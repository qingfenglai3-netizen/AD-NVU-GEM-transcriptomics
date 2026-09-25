"""Figure 1 workflow schematic assembled from the panel images."""
from pathlib import Path
import sys
R=Path(__file__).resolve().parents[1]
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch,Circle,FancyArrowPatch
from PIL import Image
A=Path(r'<ANALYSIS_ROOT>/fuxian_test/JTM_manuscript/Figure1_panel_assets')
plt.rcParams.update({'font.family':'Arial','font.size':7,'pdf.fonttype':42,'svg.fonttype':'none'})
fig=plt.figure(figsize=(10.0,5.5));ax=fig.add_axes([0,0,1,1]);ax.set(xlim=(0,1),ylim=(0,1));ax.axis('off')
navy='#05285E'
def panel(x,y,w,h,n,title,color):
 ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.009,rounding_size=0.022',fc=color,ec='#E1E6ED',lw=.5))
 ax.add_patch(Circle((x+.020,y+h-.026),.016,color=navy));ax.text(x+.020,y+h-.026,str(n),ha='center',va='center',color='white',fontweight='bold',fontsize=9)
 ax.text(x+.047,y+h-.028,title,ha='left',va='center',fontweight='bold',fontsize=10,color=navy)
def text(x,y,s,size=8,**kw):ax.text(x,y,s,ha='center',va='center',fontsize=size,**kw)
def box(x,y,w,h,s,color='#F5F8F3',size=8):
 ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=.006,rounding_size=.009',fc=color,ec='#789194',lw=.6));text(x+w/2,y+h/2,s,size,fontweight='bold')
def asset(name,rect):
 a=fig.add_axes(rect);a.imshow(Image.open(A/name));a.axis('off')
def arrow(a,b):ax.add_patch(FancyArrowPatch(a,b,arrowstyle='-|>',mutation_scale=13,lw=1.4,color='#2E3A4B'))
panel(.015,.515,.280,.464,1,'snRNA-seq datasets','#E9EEF7')
box(.040,.790,.231,.094,'GSE188545 · MTG\n4 AD / 6 CN retained',size=9)
box(.040,.663,.231,.100,'GSE237718 · temporal cortex\n29 AD / 27 CN retained',size=9)
text(.155,.582,'66 donors after QC\n389,853 nuclei',10,fontweight='bold')
panel(.015,.026,.280,.455,2,'Bulk microarray','#E9E7E7')
asset('Figure3A_bulk_volcano_full_labels.png',[.035,.113,.145,.272])
text(.232,.304,'GSE132903\nMTG\n97 AD\n98 controls',9,fontweight='bold')
box(.037,.060,.235,.043,'189 AD-up signature genes','#F7E6EB',8)
panel(.347,.026,.285,.953,3,'AD–NVU atlas and scoring','#E8F3ED')
asset('Figure2A_atlas_UMAP_transparent.png',[.358,.716,.135,.18])
asset('Figure3C_ADup_scoring_heatmap_full_labels.png',[.496,.734,.124,.146])
text(.492,.688,'Cell annotation and five-method\nAD signature scoring',9,fontweight='bold')
arrow((.492,.650),(.492,.602))
box(.369,.454,.240,.141,'NVU-focused compartments\n\nAstrocytes · Microglia\nCerebrovascular cells',size=9)
arrow((.492,.443),(.492,.390))
box(.376,.310,.226,.070,'Marker-supported subtyping',size=9)
arrow((.492,.299),(.492,.257))
box(.376,.149,.226,.096,'Donor-level MiloR\nBulk subtype-marker ssGSEA',size=9)
text(.491,.085,'Cohort and APOE sensitivity',8,fontweight='bold')
panel(.682,.515,.300,.464,4,'State and spatial context','#F6F3F7')
box(.703,.846,.259,.060,'LIANA+ communication and enrichment',size=8)
box(.703,.760,.259,.055,'scTour state ordering',size=8)
asset('Figure7C_three_spatial_maps_horizontal_full_labels_transparent.png',[.703,.600,.258,.142])
text(.832,.563,'RCTD localization · spatial LR scoring',8,fontweight='bold')
panel(.682,.026,.300,.455,5,'GEM candidate prioritization','#F8F2EC')
text(.831,.407,'Glial-ligand NicheNet → vascular GEM',8,fontweight='bold')
asset('Figure8B_GEM_virtual_KD_network_full_labels.png',[.693,.208,.134,.164])
asset('Figure9D_bulk_AUC_heatmap_full_labels.png',[.833,.207,.135,.164])
text(.756,.185,'Virtual perturbation',7,fontweight='bold');text(.900,.185,'Ancillary classification',7,fontweight='bold')
text(.832,.118,'SEA-AD vascular-cell expression\nIndependent temporal-cortex validation\nGSE5281 neuronal context',8,fontweight='bold')
text(.832,.056,'Drug-target and pathway annotation',7.5,fontweight='bold')
arrow((.306,.509),(.338,.509));arrow((.643,.754),(.673,.754))
for ext in ['png','tif','pdf','svg']:
 kw={'dpi':600,'facecolor':'white'}
 if ext=='tif':kw['pil_kwargs']={'compression':'tiff_lzw'}
 fig.savefig(R/'figures'/f'Fig1.{ext}',**kw)
plt.close(fig)
