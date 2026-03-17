import pandas as pd

# 读取文件
df_merged = pd.read_csv('merged.csv', low_memory=False)
df_sample = pd.read_csv('sample.csv')

# 时间列转datetime
df_merged['记录时间'] = pd.to_datetime(df_merged['记录时间'], errors='coerce', utc=True).dt.tz_localize(None)
df_merged['入ICU时间'] = pd.to_datetime(df_merged['入ICU时间'], errors='coerce', utc=True).dt.tz_localize(None)
df_sample['入ICU时间'] = pd.to_datetime(df_sample['入ICU时间'], errors='coerce', utc=True).dt.tz_localize(None)
df_sample['首次CRRT时间'] = pd.to_datetime(df_sample['首次CRRT时间'], errors='coerce', utc=True).dt.tz_localize(None)
df_sample['结束时间'] = pd.to_datetime(df_sample['结束时间'], errors='coerce', utc=True).dt.tz_localize(None)

# 只保留5种体征数据，并转数值
vitals_info = [
    ('T(℃)', 'T'),
    ('HR(bpm)', 'HR'),
    ('BPd(mmHg)', 'BPd'),
    ('BPs(mmHg)', 'BPs'),
    ('血糖(mmol/L)', '血糖')
]

vital_names = [v[0] for v in vitals_info]
short_names = [v[1] for v in vitals_info]

df_vitals = df_merged[df_merged['记录名称'].isin(vital_names)].copy()
df_vitals['记录值'] = pd.to_numeric(df_vitals['记录值'], errors='coerce')
df_vitals = df_vitals.dropna(subset=['记录值'])

result_rows = []

for _, sample_row in df_sample.iterrows():
    mask = (
            (df_vitals['住院号'] == sample_row['住院号']) &
            (df_vitals['年龄'] == sample_row['年龄']) &
            (df_vitals['患者姓名'] == sample_row['患者姓名']) &
            (df_vitals['性别'] == sample_row['性别']) &
            (df_vitals['入ICU时间'] == sample_row['入ICU时间'])
    )
    patient_vitals = df_vitals[mask].copy()

    if patient_vitals.empty:
        print(f"警告：患者 {sample_row['住院号']} {sample_row['患者姓名']} 无匹配体征数据，已跳过")
        continue

    crrt_time = sample_row['首次CRRT时间']
    crrt_end = sample_row['结束时间']
    icu_time = sample_row['入ICU时间']

    # 计算最大小时数
    if pd.notna(crrt_end):
        max_hours = int((crrt_end - crrt_time).total_seconds() // 3600) + 1  # 包括T0
        max_hours = min(max_hours, 73)  # 最多不超过72h
    else:
        max_hours = 73  # 默认0~72

    # 1. 入室体温：离入ICU时间最近的T
    t_rec = patient_vitals[patient_vitals['记录名称'] == 'T(℃)'].copy()
    if not t_rec.empty:
        t_rec['diff'] = (t_rec['记录时间'] - icu_time).abs()
        entry_temp = t_rec.loc[t_rec['diff'].idxmin(), '记录值']
    else:
        entry_temp = None

    # 2. CRRT前所有平均值
    before_all = patient_vitals[patient_vitals['记录时间'] < crrt_time]
    avg_before = before_all.groupby('记录名称')['记录值'].mean()

    # 3. CRRT前最后一次（最近一次）
    if not before_all.empty:
        before_last = before_all.sort_values('记录时间', ascending=False).groupby('记录名称').first()['记录值']
    else:
        before_last = pd.Series(dtype='float64')

    # 4. 小时点数据 T0 ~ T72
    hourly_data = {short: [None] * max_hours for short in short_names}

    for hour in range(max_hours):  # 0 to 72
        target_time = crrt_time + pd.Timedelta(hours=hour)
        for full_name, short in vitals_info:
            vit_rec = patient_vitals[patient_vitals['记录名称'] == full_name].copy()
            if vit_rec.empty:
                hourly_data[short].append(None)
                continue
            vit_rec['diff'] = (vit_rec['记录时间'] - target_time).abs()
            closest = vit_rec.loc[vit_rec['diff'].idxmin()]
            if closest['diff'] <= pd.Timedelta(minutes=30):
                hourly_data[short].append(closest['记录值'])
            else:
                hourly_data[short].append(None)

    # 构建结果字典
    row_dict = sample_row.to_dict()
    row_dict['入室体温'] = entry_temp

    # 添加前平均和前最后
    for full_name, short in vitals_info:
        row_dict[f"avg_before_{short}"] = avg_before.get(full_name, None)
        row_dict[f"last_before_{short}"] = before_last.get(full_name, None) if not before_last.empty else None

    # 添加小时点
    for short in short_names:
        for hour in range(max_hours):
            col_name = f"{short}{hour}"
            row_dict[col_name] = hourly_data[short][hour]

    result_rows.append(row_dict)

result_df = pd.DataFrame(result_rows)
base_cols = list(df_sample.columns) + ['入室体温']

pre_cols = []
for short in short_names:
    pre_cols += [f"avg_before_{short}", f"last_before_{short}"]

time_series_cols = []
for short in short_names:  # 按 T -> HR -> BPd -> BPs -> 血糖 的顺序
    time_series_cols += [f"{short}{h}" for h in range(73)]

final_columns = base_cols + pre_cols + time_series_cols
final_columns = [col for col in final_columns if col in result_df.columns]
result_df = result_df[final_columns]

# 输出结果
result_df.to_csv('result_grouped_by_vital.csv', index=False, encoding='utf-8-sig')

print("处理完成！结果已保存为 result_grouped_by_vital.csv")
print(f"共处理 {len(result_df)} 名患者")
print("时序列顺序示例：T0 T1 ... T72, HR0 HR1 ... HR72, BPd0 ... BPd72, BPs0 ... BPs72, 血糖0 ... 血糖72")
