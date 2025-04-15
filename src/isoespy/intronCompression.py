
def calc_populations(isoform_model):
    """ P_ranges と P_components を計算する """
    # 辞書から範囲を展開し、ソート
    indexed_ranges = [(start, end, key, i) for key, ranges in isoform_model.items() for i, (start, end) in enumerate(ranges)]
    indexed_ranges.sort()

    merged = []
    merged_indices = []
    
    for start, end, source, index in indexed_ranges:
        if (not merged) or (merged[-1][1] < start-1):
            # 新規P (現時点の最右populationとオーバーラップしないため新規)
            merged.append((start, end))
            merged_indices.append({source: [index]})
        else:
            # 既存P に追加
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
            if source in merged_indices[-1]:
                merged_indices[-1][source].append(index)
            else:
                merged_indices[-1][source] = [index]
    
    P_ranges = {i+1: merged[i] for i in range(len(merged))}
    P_components = {i+1: merged_indices[i] for i in range(len(merged))}
    for tx in isoform_model:
        for i in P_components:
            if not tx in P_components[i]:
                P_components[i][tx] = []
 
    return P_ranges, P_components


def calc_introns(p_ranges):
    sorted_ranges = sorted(p_ranges.values())
    gaps = {}
    
    gap_start = None
    gap_index = 1
    
    for i in range(len(sorted_ranges) - 1):
        gap_start = sorted_ranges[i][1] + 1
        gap_end = sorted_ranges[i+1][0] - 1
        
        if gap_start <= gap_end:
            gaps[gap_index] = (gap_start, gap_end)
            gap_index += 1
    
    return gaps


def calc_migration_distance(i_ranges, rate):
    if rate < 1:
        # 方法1: イントロン幅を 100*rate % に圧縮
        p_left = {}
        p_left[1] = 0
        for i in i_ranges:
            p_left[i+1] = 0
            # 左に詰める分だけ加算
            for j in range(1, i+1):
                p_left[i+1] += (i_ranges[j][1]-i_ranges[j][0]+1)*(1-rate)
    else:
        # 方法2: イントロン幅を rate [nt] に固定
        p_left = {}
        p_left[1] = 0
        for i in i_ranges:
            p_left[i+1] = 0
            # 左に詰める分だけ加算
            for j in range(1, i+1):
                width = i_ranges[j][1]-i_ranges[j][0]+1
                p_left[i+1] += width - rate
    return p_left


def migrate_to_left(isoform_model, p_left, p_components):
    """ isoform_model: 元の座標
        p_left:        左づめの分量
        p_components:  対象となるエクソン
    """
    moved_model = {tx: [] for tx in isoform_model}
    for tx in isoform_model:
        number_of_exons = len(isoform_model[tx])
        for i in range(number_of_exons):
            # txの第iエクソン (0, 1, 2,...)
            for j in p_components:
                # Pjに属しているか
                if i in p_components[j][tx]:
                    moved_model[tx].append((isoform_model[tx][i][0]-p_left[j], isoform_model[tx][i][1]-p_left[j]))
                    break
    return moved_model


def intronCompression(isoform_model, rate=100):
    P_ranges, P_components = calc_populations(isoform_model)
    I_ranges = calc_introns(P_ranges)
    P_left = calc_migration_distance(I_ranges, rate)
    moved_isoform_model = migrate_to_left(isoform_model, P_left, P_components)
    #MIN = float('inf')
    #MAX = -float('inf')
    #for tx in moved_isoform_model:
    #    if moved_isoform_model[tx][0][0] < MIN:
    #        MIN = moved_isoform_model[tx][0][0]
    #    if MAX < moved_isoform_model[tx][-1][1]:
    #        MAX = moved_isoform_model[tx][-1][1]
    return moved_isoform_model


if __name__ == "__main__":
    # テストケース
    isoform_model = {
        "tx1": [(30,60), (80,110), (150,180), (200,210), (250,270)],
        "tx2": [(50,90), (150,180), (230,230), (260,280)],
        "tx3": [(70,90), (150,170), (250,280)]
    }

    P_ranges, P_components = calc_populations(isoform_model)
    I_ranges = calc_introns(P_ranges)
    P_left = calc_migration_distance(I_ranges, 0.1)
    print("P_ranges", P_ranges)
    print("P_components", P_components)
    print("I_ranges", I_ranges)
    print("P_left", P_left)
    print("isoform_model", isoform_model)
    moved_isoform_model = migrate_to_left(isoform_model, P_left, P_components)
    print("moved_isoform_model", moved_isoform_model)

