"""Reproduce the effective-state example, transport counterexample and diagram."""
from pathlib import Path
from collections import defaultdict
from html import escape
import json
import math
import numpy as np

ROOT=Path(__file__).resolve().parent


def information(joint):
    px=joint.sum(1,keepdims=True)
    py=joint.sum(0,keepdims=True)
    keep=joint>0
    return float(np.sum(joint[keep]*np.log(joint[keep]/(px@py)[keep])))


def diagram():
    parts=['''<svg xmlns="http://www.w3.org/2000/svg" width="720" height="720" viewBox="0 0 720 720" role="img" aria-labelledby="title desc">
<title id="title">从有效状态到实际计算</title><desc id="desc">子树内部观测经边界状态影响外部预测；有效状态需经过统计辨认、递归复用和计算实现。</desc>
<defs><marker id="arrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse"><path d="M0 0L10 5L0 10z" fill="#6a7b91"/></marker></defs>
<rect width="720" height="720" rx="14" fill="#f8fafc"/>
<g font-family="PingFang SC, Noto Sans CJK SC, Microsoft YaHei, sans-serif" fill="#20364e">
<text x="32" y="49" font-size="28" font-weight="650">大量内部实现，少数对外预测坐标</text>
<text x="32" y="83" font-size="20" fill="#5f7085">条件分隔给出压缩的依据，统计与算法决定能否学到。</text>''']
    boxes=[(32,'内部观测 B','词形与内部实现','#e8f0fb'),(261,'边界状态 A','后验权重 q(a|b)','#e3f3ed'),(490,'外部目标 C','类别与后续词','#fff1e3')]
    for x,title,caption,color in boxes:
        parts.append(f'<rect x="{x}" y="111" width="198" height="94" rx="10" fill="{color}" stroke="#d6e0e9"/>')
        parts.append(f'<text x="{x+99}" y="149" text-anchor="middle" font-size="24" font-weight="600">{escape(title)}</text>')
        parts.append(f'<text x="{x+99}" y="181" text-anchor="middle" font-size="20">{escape(caption)}</text>')
    for x1,x2 in [(235,254),(464,483)]:
        parts.append(f'<path d="M{x1} 159H{x2}" stroke="#6a7b91" stroke-width="2" marker-end="url(#arrow)"/>')
    parts.append('''<text x="360" y="246" text-anchor="middle" font-size="26" font-family="Georgia,serif">p(c | b) = Σₐ p(c | a) p(a | b)</text>
<text x="360" y="280" text-anchor="middle" font-size="20" fill="#5f7085">给定 A 后，B 与 C 条件独立。</text>
<path d="M32 304H688" stroke="#dae3eb"/>
<text x="32" y="341" font-size="26" font-weight="600">从概率结构到学习优势，还需四个环节</text>''')
    stages=[('1  有效状态','对外关系能被少数','硬或软状态表达'),('2  统计可辨认','不同作用的差别','超过有限样本噪声'),('3  递归复用','恢复的状态继续','组织下一尺度关系'),('4  实际计算','架构与训练取得','并使用这些关系')]
    for i,(title,line1,line2) in enumerate(stages):
        x=32+344*(i%2); top=368+142*(i//2)
        parts.append(f'<rect x="{x}" y="{top}" width="312" height="119" rx="9" fill="white" stroke="#bdcddd"/>')
        for y,text,size,weight in [(top+33,title,24,600),(top+68,line1,22,400),(top+96,line2,22,400)]:
            parts.append(f'<text x="{x+156}" y="{y}" text-anchor="middle" font-size="{size}" font-weight="{weight}">{escape(text)}</text>')
        if i%2==0:parts.append(f'<path d="M{x+317} {top+57}H{x+337}" stroke="#6a7b91" stroke-width="2" marker-end="url(#arrow)"/>')
    parts.append('<text x="360" y="678" text-anchor="middle" font-size="20" fill="#5f7085">低秩不独自证明树，共同几何还需功能检验。</text></g></svg>')
    (ROOT/'figures/effective_states.svg').write_text('\n'.join(parts)+'\n')


def main():
    q=np.array([[1.,0.],[1.,0.],[0.,1.],[.5,.5]])
    r=np.array([[.8,.2],[.3,.7]])
    pb=np.array([.125,.125,.25,.5])
    m=q@r
    assert np.allclose(m,[[.8,.2],[.8,.2],[.3,.7],[.55,.45]])
    assert np.allclose(pb@q,[.5,.5])
    # Construct a full joint law independently and recover p(C|B).
    joint=pb[:,None,None]*q[:,:,None]*r[None,:,:]
    recovered=joint.sum(1)/pb[:,None]
    assert np.allclose(recovered,m)
    examples=[(prev,cur,cur+prev,.25) for prev in (-1,1) for cur in (-1,1)]
    before=np.zeros((2,2)); after=np.zeros((2,3)); residual=np.zeros((2,2))
    for prev,cur,total,prob in examples:
        a=(prev+1)//2
        before[a,(cur+1)//2]+=prob
        after[a,(total+2)//2]+=prob
        residual[a,(total-cur+1)//2]+=prob
    transport={'before':information(before),'after':information(after),'residual':information(residual)}
    assert abs(transport['before'])<1e-14
    assert abs(transport['after']-.5*math.log(2))<1e-14
    assert abs(transport['residual']-math.log(2))<1e-14
    # The positive linear-prediction residual is exactly the transported input.
    x=np.array([cur for prev,cur,total,p in examples])
    y=np.array([total for prev,cur,total,p in examples])
    a,b=np.linalg.lstsq(np.column_stack([x,np.ones(4)]),y,rcond=None)[0]
    error=np.max(np.abs(y-a*x-b-np.array([prev for prev,cur,total,p in examples])))
    assert error<1e-14
    out={'effective_state_matrix':m.tolist(),'matrix_rank':int(np.linalg.matrix_rank(m)),
         'joint_factorization_error':float(np.max(np.abs(recovered-m))),
         'transport_mutual_information_nats':transport,'transport_residual_error':float(error),
         'natural_language_exponent_arithmetic':{'TinyStories':.325/(2*.88),'WikiText':.265/(2*.94)}}
    (ROOT/'mechanism_checks.json').write_text(json.dumps(out,indent=2)+'\n')
    diagram()
    policy_path=ROOT.parents[1]/'data/hierarchical_data_figures.json'
    policy=json.loads(policy_path.read_text())
    policy['figures/effective_states']={'width':900,'kind':'original-mechanism-diagram','intrinsic_width':720,'intrinsic_height':720}
    policy_path.write_text(json.dumps(policy,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':main()
