"""Three-page, bilingual TagForge figures: self-contained HTML -> vector PDF.

Page 1: generic workflow and extensible application scenarios.
Page 2: correction decisions, UMI graph rules and implementation-backed strengths.
Page 3: PI-seq molecule assignment and PB–CB cell mapping.
Run with Python + matplotlib + fontTools and Google Chrome.
English: Arial. Chinese: embedded Droid Sans Fallback. Counts are illustrative.
"""
from pathlib import Path
import base64
from io import BytesIO, StringIO
import subprocess
import tempfile
from math import hypot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.path import Path as DrawingPath
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Rectangle
from fontTools import subset
from fontTools.ttLib import TTFont

OUT = Path(__file__).resolve().parent
FONT_DIR = Path('/data/project/yanhuixu/yfu/01_source/fonts')
REGULAR, BOLD, CJK = [FONT_DIR / n for n in ('Arial.ttf', 'Arial Bold.ttf', 'Arial Unicode.ttf')]
# Arial Unicode is only used to measure the old layout. Browser rendering uses
# Arial for Latin glyphs and the dedicated CJK font for Chinese glyphs.
CJK_RENDER = Path('/usr/share/fonts/truetype/droid/DroidSansFallbackFull.ttf')
for path in (REGULAR, BOLD, CJK):
    font_manager.fontManager.addfont(str(path))
plt.rcParams.update({'pdf.fonttype': 42, 'ps.fonttype': 42, 'pdf.compression': 9,
                     'svg.fonttype': 'none'})
NAVY, MUTED, LINE = '#142D46', '#5B6D7C', '#D7E0E7'
BLUE, TEAL, PURPLE, RED = '#2876B6', '#12877D', '#7860AE', '#BB5962'
PALE_BLUE, PALE_TEAL, PALE_PURPLE = '#EDF5FC', '#EAF6F3', '#F3EFF9'


class Page:
    def __init__(self, language):
        self.zh = language == 'zh-CN'
        self.fig = plt.figure(figsize=(16, 10), facecolor='white')
        self.ax = self.fig.add_axes([0, 0, 1, 1])
        self.ax.set(xlim=(0, 1440), ylim=(0, 900))
        self.ax.axis('off')
        self.labels = []
        self.width_checks = []

    def tr(self, en, cn):
        return cn if self.zh else en

    def txt(self, x, y, label, size=12, color=NAVY, bold=False, ha='left', max_width=None):
        has_cjk = any('\u3400' <= c <= '\u9fff' for c in label)
        font = CJK if has_cjk else BOLD if bold else REGULAR
        obj = self.ax.text(x, y, label, fontsize=size, color=color,
                           fontproperties=font_manager.FontProperties(fname=str(font)),
                           ha=ha, va='center', linespacing=1.5, zorder=6)
        self.labels.append(obj)
        if bold:
            obj.set_gid(f'bold-label-{len(self.labels)}')
        if max_width:
            self.width_checks.append((obj, max_width))
        return obj

    def box(self, x, y, w, h, face='white', edge=LINE, radius=12, lw=1):
        p = FancyBboxPatch((x, y), w, h, boxstyle=f'round,pad=0,rounding_size={radius}',
                          facecolor=face, edgecolor=edge, linewidth=lw, zorder=1)
        self.ax.add_patch(p)

    def line(self, points, color=LINE, lw=1.4, zorder=2):
        self.ax.plot(*zip(*points), color=color, lw=lw, zorder=zorder,
                     solid_capstyle='round', solid_joinstyle='round')

    def arrow(self, start, end, color=MUTED, lw=1.4, scale=9):
        self.ax.add_patch(FancyArrowPatch(start, end, arrowstyle='-|>',
                         mutation_scale=scale, color=color, linewidth=lw,
                         shrinkA=0, shrinkB=0, connectionstyle='arc3,rad=0',
                         capstyle='round', joinstyle='round', zorder=3))

    def route(self, points, color=MUTED, lw=1.4, scale=9, radius=8):
        """Draw one continuous rounded path, with its arrowhead at the endpoint."""
        vertices=[points[0]]
        codes=[DrawingPath.MOVETO]
        for before,corner,after in zip(points,points[1:],points[2:]):
            incoming=hypot(corner[0]-before[0],corner[1]-before[1])
            outgoing=hypot(after[0]-corner[0],after[1]-corner[1])
            r=min(radius,incoming/2,outgoing/2)
            entry=tuple(corner[j]+(before[j]-corner[j])*r/incoming for j in (0,1))
            leave=tuple(corner[j]+(after[j]-corner[j])*r/outgoing for j in (0,1))
            vertices.extend([entry,corner,leave])
            codes.extend([DrawingPath.LINETO,DrawingPath.CURVE3,DrawingPath.CURVE3])
        vertices.append(points[-1])
        codes.append(DrawingPath.LINETO)
        self.ax.add_patch(FancyArrowPatch(path=DrawingPath(vertices,codes),
                         arrowstyle='-|>', mutation_scale=scale, color=color,
                         linewidth=lw, shrinkA=0, shrinkB=0,
                         capstyle='round', joinstyle='round', zorder=3))

    def pill(self, x, y, w, label, face, color, size=11, h=27):
        self.box(x, y, w, h, face=face, edge=face, radius=6)
        self.txt(x+w/2, y+h/2, label, size, color, ha='center', max_width=w-8)

    def bases(self, x, y, sequence, highlight=(), size=14, step=17):
        for i, base in enumerate(sequence):
            if i in highlight:
                self.box(x+i*step-7, y-12, 15, 24, '#FBECEE', '#FBECEE', 3)
            self.txt(x+i*step, y, base, size, RED if i in highlight else NAVY, ha='center')

    def node(self, x, y, r, color, umi, count):
        self.ax.add_patch(Circle((x, y), r, facecolor=color, edgecolor='white', lw=1.2, zorder=4))
        self.txt(x, y+6, umi, 10.5, 'white', True, 'center')
        self.txt(x, y-11, f'{count} {"read" if count == 1 else "reads"}', 9, 'white', ha='center')

    def header(self, title, subtitle, page, extension=False):
        self.box(48, 826, 7, 45, TEAL, TEAL, 2)
        self.txt(72, 849, 'TagForge', 34, bold=True)
        self.txt(1392, 850, self.tr('PI-seq  /  EXPERIMENT-SPECIFIC EXTENSIONS',
                                  'PI-seq  /  实验定制扩展') if extension else
                 self.tr('GENERAL CORE  /  BARCODE-TO-MOLECULE ANALYSIS',
                         '通用核心  /  从条码证据到分子计数'), 11, MUTED, ha='right')
        self.txt(48, 793, title, 23, bold=True, max_width=1340)
        self.txt(48, 758, subtitle, 13, MUTED, max_width=1340)
        self.txt(1392, 28, page, 9, MUTED, ha='right')

    def panel(self, x, y, w, h, number, title, color, face='white'):
        self.box(x, y, w, h, face, face if face != 'white' else LINE)
        self.pill(x+17, y+h-42, 39, number, color, 'white', 12)
        self.txt(x+68, y+h-29, title, 16, bold=True, max_width=w-85)
        self.line([(x+18, y+h-53), (x+w-18, y+h-53)], lw=.8)

    def finish(self, pdf):
        self.fig.canvas.draw()
        renderer = self.fig.canvas.get_renderer()
        canvas = self.fig.bbox
        for label in self.labels:
            bb = label.get_window_extent(renderer)
            if bb.x0 < 0 or bb.y0 < 0 or bb.x1 > canvas.x1 or bb.y1 > canvas.y1:
                raise RuntimeError(f'Label outside canvas: {label.get_text()}')
        for label, width in self.width_checks:
            bb = label.get_window_extent(renderer)
            actual = bb.width / canvas.width * 1440
            if actual > width:
                raise RuntimeError(f'Label too wide ({actual:.1f} > {width}): {label.get_text()}')
        buffer = StringIO()
        self.fig.savefig(buffer, format='svg', metadata={'Date': None})
        svg = buffer.getvalue()
        svg = svg[svg.index('<svg'):]
        pdf.append((svg, '\n'.join(label.get_text() for label in self.labels)))
        plt.close(self.fig)


def overview(language, pdf):
    p = Page(language)
    t = p.tr
    p.header(t('From complex barcode sequencing to molecular quantification', '从复杂条码测序数据到可追踪的分子证据与定量结果'),
             t('Parse structure → correct barcodes → count molecules',
               '结构解析 → 条码校正 → 分子计数'), '1 / 3')
    p.panel(48, 480, 416, 246, '01', t('Parse library structure', '解析文库结构'), BLUE)
    p.panel(500, 480, 416, 246, '02', t('Correct each segment', '逐片段校正条码'), PURPLE)
    p.panel(952, 480, 440, 246, '03', t('UMI clustering & counting', 'UMI 聚类与分子计数'), TEAL)
    p.arrow((470, 606), (494, 606), BLUE)
    p.arrow((922, 606), (946, 606), PURPLE)
    p.txt(69, 651, t('B1: identity (e.g. cell / sample); B2: feature', 'B1：身份条码（如细胞 / 样本）；B2：特征条码'), 10.1, MUTED, max_width=374)
    p.txt(69, 614, 'R1', 12, bold=True)
    p.line([(108,614), (439,614)], lw=1.2, zorder=0)
    for x,w,label,face,color in [(108,48,'L', '#EEF1F5',MUTED), (163,102,'B1a',PALE_BLUE,BLUE),
                                (272,48,'L','#EEF1F5',MUTED), (327,102,'B1b',PALE_BLUE,BLUE)]:
        p.box(x,601,w,26,face,color,4,.7)
        p.txt(x+w/2,614,label,11,color,ha='center')
    p.txt(69, 568, 'R2', 12, bold=True)
    p.line([(108,568), (432,568)], lw=1.2, zorder=0)
    for x,w,label,face,color in [(108,88,'U1',PALE_TEAL,TEAL),(206,128,'B2',PALE_PURPLE,PURPLE),
                                (344,88,'U2',PALE_TEAL,TEAL)]:
        p.box(x,555,w,26,face,color,4,.7)
        p.txt(x+w/2,568,label,11,color,ha='center')
    p.txt(69,524,t('Linker localization; fixed fallback if configured', 'Linker 定位 / 固定位置（按配置回退）'),11.3)
    p.txt(521,651,t('Per-segment whitelist correction', '逐片段白名单校正'),11.5,MUTED)
    p.txt(521,615,t('Mismatch','碱基错配'),12)
    p.bases(635,615,'ACGA',{3})
    p.arrow((704,615),(741,615),PURPLE,lw=1.2,scale=8.5)
    p.bases(754,615,'ACGT')
    p.txt(521,574,t('Shift','位置偏移'),12)
    p.bases(635,574,'TACGTG',{0,5},13,14)
    p.arrow((720,574),(746,574),PURPLE,lw=1.2,scale=8.5)
    p.bases(760,574,'ACGT',size=13,step=15)
    p.txt(521,524,t('Correct segments, then assemble barcodes', '校正后拼接条码身份'),11.3,PURPLE)
    p.txt(973,651,t('Group by identity barcode × feature', '按身份条码 × 特征分组'),11.3,MUTED)
    p.box(983,559,123,62,'#F7F9FB','#F7F9FB',9)
    p.txt(1044,590,'12',27,bold=True,ha='center')
    p.txt(1044,540,'reads',12,MUTED,ha='center')
    p.arrow((1121,590),(1190,590),TEAL,scale=10)
    p.box(1207,559,151,62,PALE_TEAL,PALE_TEAL,9)
    p.txt(1282,590,'2',27,TEAL,bold=True,ha='center')
    p.txt(1282,540,t('molecules','分子'),12,TEAL,ha='center')
    # Only the generic output is prominent. The diagnostic path starts from it.
    p.route([(1172,474),(1172,452),(371,452),(371,430)],TEAL,radius=9)
    p.panel(48,213,646,211,'04',t('Molecule evidence → sparse counts','分子证据 → 稀疏计数'),BLUE,PALE_BLUE)
    p.panel(742,213,650,211,'05',t('Downsampling diagnostics','降采样诊断'),PURPLE,PALE_PURPLE)
    p.arrow((700,307),(736,307),BLUE)
    # Molecular detail table: all columns are retained by the generic engine.
    xs=[75,145,238,340]
    for xx, label in zip(xs,['B1','Feature','UMI','Reads']):
        p.txt(xx,346,label,9.5,MUTED)
    for yy, row in [(317,['B01','F1','AAAC','9']), (287,['B01','F1','GGGC','3'])]:
        for xx,label in zip(xs,row):
            p.txt(xx,yy,label,11)
    p.line([(69,332),(423,332)],'white',1)
    p.arrow((432,304),(474,304),BLUE,lw=1.2,scale=8.5)
    p.box(482,276,64,57,'white','white',7)
    p.txt(514,304,'2',24,BLUE,True,'center')
    p.txt(568,317,'B01 × F1',11,BLUE)
    p.txt(568,288,t('molecules','个分子'),11,MUTED)
    p.txt(69,244,'Parquet · MTX · TSV',11,BLUE)
    p.txt(763,346,t('Downsample molecule-supporting reads', '对分子支持 reads 降采样'),11.2,MUTED)
    p.txt(763,305,t('Sequencing saturation','测序饱和度'),13,BLUE)
    p.txt(763,276,t('Read redundancy','reads 的重复程度'),11,MUTED)
    p.txt(1090,305,t('Non-singleton fraction','非单例分子比例'),12.5,RED,max_width=281)
    p.txt(1090,276,t('Fraction with ≥2 reads','多条 reads 支持的分子占比'),11,MUTED,max_width=281)
    # Applications share the same barcode x feature x UMI input model.
    # Scale/recovery mechanisms are detailed on page 2.
    p.txt(48,181,t('APPLICATIONS & EXTENSIONS', '应用与扩展场景'),14,bold=True)
    p.txt(365,181,t('Identity barcode × feature barcode × UMI',
                    '身份条码 × 特征条码 × UMI'),10.6,MUTED,max_width=760)
    p.pill(1180,168,212,t('CONFIGURE · EXTEND', '结构可配置 · 场景可扩展'),
           PALE_TEAL,TEAL,10.3)
    scenarios = [
        (t('Antibody & combinatorial tags','抗体标签与组合条码'),
         t('Identity barcode × antibody tag','身份条码 × 抗体标签'),BLUE,PALE_BLUE),
        ('Sample Tags',
         t('Cell × sample tag','细胞 × 样本标签'),PURPLE,'white'),
        ('CRISPR Guide Capture',
         t('Cell × sgRNA','细胞 × sgRNA'),BLUE,'white'),
        ('Antigen Capture',
         t('Cell × antigen tag','细胞 × 抗原标签'),PURPLE,'white'),
        (t('Promoter / reporter barcodes','启动子 / 报告基因条码'),
         t('Identity barcode × design tag','身份条码 × 设计标签'),TEAL,PALE_TEAL),
    ]
    for i,(title,output,color,face) in enumerate(scenarios):
        x=48+i*273
        p.box(x,76,252,77,face,LINE if face=='white' else face,9,.8)
        p.box(x+12,94,3,40,color,color,1)
        p.txt(x+25,126,title,11.1,bold=True,max_width=215)
        p.txt(x+25,96,output,10.5,color,max_width=215)
    p.txt(48,28,t('L: linker · U: UMI. Illustrative layout and counts; extensions require validation.',
                  'L：linker · U：UMI；结构与计数为示意，扩展场景需验证。'),9.2,MUTED)
    p.finish(pdf)


def mechanisms(language,pdf):
    p=Page(language)
    t=p.tr
    p.header(t('Core mechanisms: barcode identity and UMI-defined molecules',
               '核心机制：从条码身份判定到 UMI 分子归并'),
             t('Whitelist correction resolves identity; UMI clustering collapses reads.',
               '白名单校正条码，UMI 聚类归并 reads。'), '2 / 3')
    p.panel(48,395,646,331,'A',t('Barcode identity correction','条码身份校正'),PURPLE)
    p.panel(742,395,650,331,'B',t('UMI molecule grouping','UMI 分子归并'),TEAL)

    # A. Keep the decision and ranking order; the full equation is in the notes.
    p.pill(69,600,202,t('Original-window exact hit?','原位精确命中？'),PALE_PURPLE,PURPLE,11.5,h=34)
    p.arrow((276,617),(335,617),PURPLE)
    p.txt(303,635,t('YES','命中'),9.5,PURPLE,ha='center')
    p.pill(340,600,332,t('Retain this barcode identity','直接保留该条码身份'),PALE_PURPLE,PURPLE,12,h=34)
    p.arrow((170,597),(170,572),MUTED)
    p.txt(194,585,t('No exact hit','未命中'),10,MUTED)
    p.pill(69,535,603,t('Search allowed shifted windows and whitelist mismatches',
                        '在允许范围内搜索偏移窗口与白名单错配候选'),
           '#F7F9FB',NAVY,11.5,h=34)
    p.arrow((370,532),(370,514),PURPLE,scale=8)
    p.box(69,479,603,32,'#F7F9FB','#F7F9FB',6)
    p.txt(83,495,t('Rank by','依次比较'),9.8,MUTED)
    for i,label in enumerate([
        t('1 Edit types','1 修改类型数'),t('2 Total change','2 总修改幅度'),
        t('3 Mismatches','3 错配数'),t('4 Shift size','4 偏移幅度')]):
        x=158+i*128
        p.pill(x,480,111,label,PALE_PURPLE,PURPLE,10,h=30)
    p.line([(370,477),(370,467)],PURPLE,1.3)
    p.route([(370,467),(213,467),(213,451)],PURPLE,lw=1.3,scale=8,radius=4)
    p.route([(370,467),(523,467),(523,451)],RED,lw=1.3,scale=8,radius=4)
    p.pill(69,416,288,t('Unique best identity → retain','唯一最优身份 → 保留'),
           PALE_PURPLE,PURPLE,11,h=32)
    p.pill(375,416,297,t('No candidate / tied identities → exclude','无候选 / 身份歧义 → 排除'),
           '#FBECEE',RED,10.4,h=32)

    # B. Two independent groups isolate the effect of read-count evidence.
    p.pill(763,594,608,t('Scope: corrected identity barcode × feature',
                        '聚类范围：校正后的身份条码 × 特征'),PALE_TEAL,TEAL,12,h=34)
    p.pill(763,548,608,'Hamming ≤ 1; n(parent) ≥ 2 × n(child) − 1',PALE_TEAL,TEAL,12,h=30)
    for x,child,merge,key in [(763,1,True,'B01 × F1'),(1075,7,False,'B01 × F2')]:
        p.box(x,420,296,112,'#F7F9FB','#F7F9FB',8)
        p.txt(x+16,515,key,10,MUTED)
        p.node(x+48,477,25,TEAL,'AAAC',8)
        p.node(x+245,477,25,'#7DBBB0','AAAT',child)
        p.txt(x+147,497,'1 bp',9.5,MUTED,ha='center')
        if merge:
            # The UMI-tools directional edge runs from parent to child.
            p.arrow((x+77,477),(x+216,477),TEAL)
        else:
            p.line([(x+80,477),(x+127,477)],'#BCC9D5',1.2)
            p.line([(x+167,477),(x+214,477)],'#BCC9D5',1.2)
            p.txt(x+147,477,'×',19,RED,ha='center')
        p.txt(x+148,438,t('9 reads → 1 molecule','9 条 reads → 1 个分子') if merge else
                            t('15 reads → 2 molecules','15 条 reads → 2 个分子'),
              11,TEAL if merge else BLUE,ha='center',max_width=272)

    # Engineering capabilities: each claim describes an implemented mechanism.
    p.txt(48,366,t('ENGINEERING ADVANTAGES','工程化优势'),14,bold=True)
    capabilities=[
        (t('Scale and optimize computation','规模化处理与计算优化'),
         [t('Block reads · Parquet · DuckDB · sparse output', '分块读取 · Parquet · DuckDB · 稀疏输出'),
          t('Exact matching · neighbor lookup · correction cache', '精确匹配 · 邻域查找 · 校正缓存')],BLUE,PALE_BLUE),
        (t('Parallelism across independent groups','按独立分组并行'),
         [t('Bucket complete barcode–feature pairs; preserve UMI clustering boundaries.',
            '完整条码–特征组合分桶，保留 UMI 聚类边界。')],PURPLE,PALE_PURPLE),
        (t('Recover and validate results','断点恢复与结果校验'),
         [t('Resume committed extraction / correction batches.', '提取 / 校正批次恢复'),
          t('Validate inputs, config and upstream state before reuse.', '复用前核验输入、配置与上游状态')],TEAL,PALE_TEAL),
        (t('Batch execution and progress tracking','批量运行与过程追踪'),
         [t('Slurm · live progress · segment QC', 'Slurm · 实时进度 · 分段质控'),
          t('Sample and batch HTML / Excel reports', '样本与批次 HTML / Excel 报告')],BLUE,PALE_BLUE),
    ]
    for i,(title,body,color,face) in enumerate(capabilities):
        x=48+(i%2)*694
        y=221-(i//2)*136
        p.box(x,y,650 if i%2 else 646,119,face,face,10)
        p.box(x+18,y+22,3,77,color,color,1)
        p.txt(x+34,y+88,title,14,bold=True,max_width=590)
        ys=(y+55,y+29) if len(body)==2 else (y+43,)
        for yy,label in zip(ys,body):
            p.txt(x+34,yy,label,11.6,MUTED,max_width=590)
    p.txt(48,28,t('UMI-tools directional, default distance 1; illustrative counts.',
                  'UMI-tools directional（默认距离 1）；计数为示意。'),9.2,MUTED)
    p.finish(pdf)


def experimental_extensions(language, pdf):
    p = Page(language)
    t = p.tr
    p.header(t('Two dedicated libraries: molecular identity and cell linkage',
               '两类定制文库：分子归属与细胞映射'),
             t('PB: shared identity barcode; FB: antibody feature barcode; CB: cell barcode.',
               'PB：两类文库共享的身份条码；FB：抗体特征条码；CB：细胞条码。'),
             '3 / 3', extension=True)
    p.pill(48, 687, 1344,
           t('Shared core: structure parsing → barcode handling → UMI clustering → supporting reads',
             '复用通用核心：结构解析 → 条码处理 → UMI 聚类 → 保留支持 reads'),
           '#F0F4F7', NAVY, 12.2, h=39)
    p.panel(48, 214, 646, 452, 'A', t('PI-seq: assign molecules to PBs', 'PI-seq：确定分子的 PB 归属'), PURPLE)
    p.panel(742, 214, 650, 452, 'B', t('PB–CB: map PBs to cells', 'PB–CB：将 PB 映射到细胞'), TEAL)
    p.txt(69, 582, t('Which PB carries this antibody-tagged molecule?',
                     '同一抗体特征与 UMI，应该归属于哪个 PB？'), 12.5, PURPLE, max_width=601)
    p.txt(763, 582, t('Which cell barcode is supported by this PB?',
                      '一个 PB 应该映射到哪个细胞条码 CB？'), 12.5, TEAL, max_width=607)
    # The PI-seq layout is the development-fixture example, not a universal layout.
    p.txt(69, 541, 'R1', 11, bold=True)
    for i in range(5):
        x=110+i*111
        p.pill(x, 527, 77, f'PB{i+1}', PALE_PURPLE, PURPLE, 11, h=28)
        if i<4:
            p.txt(x+94,541,'L',10,MUTED,ha='center')
    p.txt(69, 500, 'R2', 11, bold=True)
    for x,w,label,face,color in [(110,130,'U1 (8 nt)',PALE_BLUE,BLUE),
                                (254,182,'FB (12 nt)',PALE_PURPLE,PURPLE),
                                (450,130,'U2 (8 nt)',PALE_BLUE,BLUE)]:
        p.pill(x,486,w,label,face,color,11,h=28)
    p.pill(763,527,155,'PB',PALE_PURPLE,PURPLE,11,h=28)
    p.txt(932,541,'×',12,MUTED,ha='center')
    p.pill(946,527,197,t('Observed CB','实测 CB 序列'),PALE_TEAL,TEAL,11,h=28)
    p.txt(1157,541,'×',12,MUTED,ha='center')
    p.pill(1171,527,200,'UMI',PALE_BLUE,BLUE,11,h=28)
    p.txt(763,500,t('No CB whitelist; optional conservative observed-CB correction.',
                     'CB 无需白名单；可按观测支持进行保守纠错。'),11.2,MUTED,max_width=606)
    p.pill(69,451,603,t('Within each FB × corrected UMI: compare PB read support',
                        '按 FB × 校正 UMI 分组，比较各 PB 的 read 支持'),PALE_PURPLE,PURPLE,11.5,h=30)
    p.pill(763,451,608,t('For each PB: aggregate CB read support after CB handling',
                         '按 PB 汇总 CB 处理后的 read 支持，选择优势 CB'),PALE_TEAL,TEAL,11.5,h=30)
    # Identical 9:1 vs 8:2 examples make the strict threshold explicit.
    for origin,color in [(69,PURPLE),(763,TEAL)]:
        is_pi=origin==69
        for j in range(2):
            x=origin+j*310
            p.box(x,321,293 if is_pi else 298,111,'#F7F9FB','#F7F9FB',8)
            p.txt(x+15,412,('FB1 × UMI-a' if j==0 else 'FB2 × UMI-b') if is_pi else
                  ('PB01' if j==0 else 'PB02'),10,MUTED)
            p.txt(x+15,380,('PB01 : PB02' if is_pi else 'CB01 : CB02')+
                  (' = 9 : 1' if j==0 else ' = 8 : 2'),11.5,max_width=270)
            label=(t('90% → retain PB01','90% → 保留 PB01') if is_pi else
                   t('90% → map to CB01','90% → 映射到 CB01')) if j==0 else \
                  (t('80% → exclude group','80% → 排除整组') if is_pi else
                   t('80% → no assignment','80% → 不分配 CB'))
            p.txt(x+15,346,label,11.5,color if j==0 else RED,max_width=270)
    p.txt(69,291,t('Single-PB groups are retained directly.',
                    '仅有一个 PB 的组直接保留。'),10.8,MUTED,max_width=600)
    p.txt(763,291,t('Read support selects the CB; UMI counts remain available.',
                     'read 支持决定归属，同时保留 UMI 分子计数。'),10.8,MUTED,max_width=606)
    p.txt(69,261,t('Output: PI-filtered PB × FB molecule details',
                    '输出：PI 过滤后的 PB × FB 分子明细'),12,PURPLE,bold=True,max_width=600)
    p.txt(69,235,'molecule_detail.rmMP.parquet',10.5,MUTED)
    p.txt(763,261,t('Output: PB → CB map and PB counts per cell',
                     '输出：PB → CB 映射及每个细胞的 PB 数'),12,TEAL,bold=True,max_width=606)
    p.txt(763,235,'pb_cb_map.tsv.gz · cb_pb_counts.tsv.gz',10.5,MUTED)
    p.txt(48,183,t('LINK THE TWO RESULTS THROUGH SHARED PBs',
                    '通过共享 PB 连接两类结果'),13,bold=True)
    for x,w,label,face,color in [
        (48,390,t('PI-seq: PB × antibody feature','PI-seq：PB × 抗体特征'),PALE_PURPLE,PURPLE),
        (504,330,t('PB–CB: PB → cell','PB–CB：PB → 细胞'),PALE_TEAL,TEAL),
        (900,492,t('Downstream: cell-associated feature evidence','下游：关联到细胞的抗体特征证据'),PALE_BLUE,BLUE)]:
        p.pill(x,122,w,label,face,color,12,h=39)
    p.txt(471,141,'+',22,MUTED,ha='center')
    p.arrow((842,141),(892,141),TEAL)
    p.txt(48,91,t('The two outputs can be joined on PB; downstream integration is performed separately.',
                   '两类输出可按 PB 关联，用于后续整合分析。'),10.8,MUTED,max_width=1344)
    p.txt(48,28,t('Examples only. Dominant reads / total reads > threshold (default 0.8); layout is configurable.',
                  '计数为示意；优势 reads / 总 reads > 阈值（默认 0.8）；具体片段结构按文库配置。'),9.2,MUTED)
    p.finish(pdf)


def embedded_font(path, characters):
    font = TTFont(path)
    options = subset.Options()
    options.drop_tables += ['FFTM']
    sub = subset.Subsetter(options=options)
    sub.populate(text=characters)
    sub.subset(font)
    font.flavor = 'woff'
    buffer = BytesIO()
    font.save(buffer)
    return base64.b64encode(buffer.getvalue()).decode('ascii')


def export_html_pdf(language, pages):
    characters = ''.join(text for _, text in pages)
    covered = set()
    for path in (REGULAR, BOLD, CJK_RENDER):
        covered.update(TTFont(path).getBestCmap())
    missing = {c for c in characters if not c.isspace() and ord(c) not in covered}
    if missing:
        raise RuntimeError(f'Font coverage missing: {missing}')
    fonts = []
    for name, weight, path in [('TF Arial', 400, REGULAR), ('TF Arial', 700, BOLD),
                               ('TF CJK', 400, CJK_RENDER)]:
        fonts.append(f'''@font-face {{font-family:"{name}"; font-style:normal;
            font-weight:{weight}; font-display:block;
            src:url(data:font/woff;base64,{embedded_font(path, characters)}) format("woff");}}''')
    html_path = OUT / f'TagForge_core_principle.{language}.html'
    pdf_path = html_path.with_suffix('.pdf')
    sections = '\n'.join(f'<section class="sheet" aria-label="Page {i+1}">{svg}</section>'
                         for i, (svg, _) in enumerate(pages))
    html = f'''<!doctype html>
<html lang="{language}"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>TagForge — Core methods and experimental extensions</title>
<style>{''.join(fonts)}
* {{box-sizing:border-box;}}
html,body {{margin:0;padding:0;background:#edf1f5;}}
.sheet {{width:min(1536px,calc(100vw - 32px));aspect-ratio:8/5;margin:24px auto;
         background:white;box-shadow:0 3px 22px #142d4615;}}
.sheet svg {{width:100%;height:100%;display:block;}}
svg text {{font-family:"TF Arial","TF CJK",sans-serif !important;
           font-synthesis:weight; text-rendering:optimizeLegibility;}}
g[id^="bold-label-"] text {{font-weight:700 !important;}}
@page {{size:16in 10in;margin:0;}}
@media print {{
 html,body {{background:white;width:16in;}}
 .sheet {{width:16in;height:10in;aspect-ratio:auto;margin:0;box-shadow:none;
          break-inside:avoid;break-after:page;}}
 .sheet:last-of-type {{break-after:auto;}}
 * {{-webkit-print-color-adjust:exact;print-color-adjust:exact;}}
}}
</style></head><body>
{sections}
<script>
document.fonts.ready.then(() => {{
  const sheets = [...document.querySelectorAll('.sheet')];
  const texts = [...document.querySelectorAll('svg text')];
  const escaped = texts.filter(node => {{
    const bounds = node.getBoundingClientRect();
    const sheet = node.closest('.sheet').getBoundingClientRect();
    return bounds.left < sheet.left-1 || bounds.right > sheet.right+1 ||
           bounds.top < sheet.top-1 || bounds.bottom > sheet.bottom+1;
  }}).map(node => node.textContent);
  document.body.dataset.layoutCheck = JSON.stringify({{
    pages:sheets.length,texts:texts.length,escaped,
    fonts:[...document.fonts].map(f=>({{family:f.family,weight:f.weight,status:f.status}}))
  }});
}});
</script></body></html>'''
    html_path.write_text('\n'.join(line.rstrip() for line in html.splitlines())+'\n')
    with tempfile.TemporaryDirectory(prefix='tagforge-figure-chrome-') as profile:
        command = ['google-chrome', '--headless=new', '--no-sandbox', '--disable-gpu',
                   '--disable-dev-shm-usage', '--disable-background-networking',
                   '--no-first-run', '--no-pdf-header-footer', '--allow-file-access-from-files',
                   f'--user-data-dir={profile}', '--virtual-time-budget=5000',
                   f'--print-to-pdf={pdf_path}', html_path.as_uri()]
        result = subprocess.run(command, capture_output=True, text=True, timeout=60)
        if result.returncode or not pdf_path.exists():
            raise RuntimeError(result.stderr[-3000:])
    print(html_path)
    print(pdf_path)


if __name__ == '__main__':
    for language in ('zh-CN','en'):
        pages=[]
        overview(language,pages)
        mechanisms(language,pages)
        experimental_extensions(language,pages)
        export_html_pdf(language,pages)
